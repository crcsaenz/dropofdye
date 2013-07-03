#include "stdafx.h"
#include "integrators/volumephoton.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"

/* Photon structs from photonmap.cpp */
struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};

struct PhotonProcess {
    // PhotonProcess Public Methods defined in photonmap.cpp
    PhotonProcess(uint32_t mp, ClosePhoton *buf);
    void operator()(const Point &p, const Photon &photon, float dist2,
         float &maxDistSquared);
    ClosePhoton *photons;
    uint32_t nLookup, nFound;
};


struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, float md2 = INFINITY)
        : photon(p), distanceSquared(md2) { }
    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
            (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }
    const Photon *photon;
    float distanceSquared;
};

/* Copied from photonmap.cpp */
class VolumePhotonShootingTask : public Task {
public:
    VolumePhotonShootingTask(int tn, float ti, Mutex &m, VolumePhotonIntegrator *in,
        ProgressReporter &prog, bool &at, int &ndp,
        vector<Photon> &volume, vector<Spectrum> &rpR, vector<Spectrum> &rpT,
		uint32_t &ns, Distribution1D *distrib, const Scene *sc, const Renderer *sr,
		float sz)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
      abortTasks(at), nDirectPaths(ndp),
      volumePhotons(volume), rpReflectances(rpR), rpTransmittances(rpT),
      nshot(ns), lightDistribution(distrib), scene(sc), renderer(sr), stepSize(sz) { }
    void Run();

    int taskNum;
    float time;
    Mutex &mutex;
    VolumePhotonIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    int &nDirectPaths;
    vector<Photon> &volumePhotons;
    vector<Spectrum> &rpReflectances, &rpTransmittances;
    uint32_t &nshot;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
	float stepSize;
};

/* Copied from photonmap.cpp */
inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot)
{
    return (found < needed && (found == 0 || found < shot / 1024));
}

/* Copied from photonmap.cpp */
inline float kernel(const Photon *photon, const Point &p, float maxDist2)
{
    float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * INV_PI * s * s;
}

/* Modified from photonmap.cpp */
Spectrum LVolumePhoton(KdTree<Photon> *map, int nPaths, int nLookup,
		ClosePhoton *lookupBuf, VolumeRegion *vr, const Point &isectPt, 
		const Vector &wo, float maxDist2)
{
	// From LPhoton
	Spectrum L(0.);
	PhotonProcess proc(nLookup, lookupBuf);
	map->Lookup(isectPt, proc, maxDist2);
	if (!proc.nFound) return Spectrum(0.);

	ClosePhoton *photons = proc.photons;
	int nFound = proc.nFound;
	for (int i = 0; i < nFound; ++i) {
		const Photon *p = photons[i].photon;
		L += vr->p(isectPt, p->wi, wo, 0.1f) * p->alpha;
	}
	
	// From EPhoton, with spherical volume instead of circular area
	return L / (nPaths * 4. / 3. * M_PI * pow(maxDist2, 1.5));
}


/* Modified from photonmap.cpp */
VolumePhotonIntegrator::VolumePhotonIntegrator(int nvol, int nl, 
		int mdepth, int mphodepth, float mdist, bool fg, float ga,
        float ss, float cs)
{
	nVolumePhotonsWanted = nvol;
	nLookup = nl;
	maxSpecularDepth = mdepth;
	maxPhotonDepth = mphodepth;
	maxDistSquared = mdist * mdist;
	finalGather = fg;
	cosGatherAngle = cos(Radians(ga));
	nVolumePaths = 0;
	volumeMap = NULL;
    
    stepSize = ss;
	causticScale = cs;
}

VolumePhotonIntegrator::~VolumePhotonIntegrator()
{
	delete volumeMap;
}

/* Copied from single.cpp */
void VolumePhotonIntegrator::RequestSamples(Sampler *sampler,
		Sample *sample, const Scene *scene)
{
	tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}

/* Copied from single.cpp with modification from photonmap.cpp */
Spectrum VolumePhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena) const
{
	VolumeRegion *vr = scene->volumeRegion;
    float t0, t1;
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }
    // Do single scattering volume integration in _vr_
    Spectrum Lv(0.);

    // Prepare for volume integration stepping
    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;

    // Compute sample patterns for single scattering samples
    float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        // Advance to sample at _t0_ and update _T_
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);

        // Possibly terminate ray marching if transmittance is small
        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        // Compute single-scattering source term at _p_
        Lv += Tr * vr->Lve(p, w, ray.time);
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && scene->lights.size() > 0) {
            int nLights = scene->lights.size();
            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
                         nLights-1);
            Light *light = scene->lights[ln];
            // Add contribution of _light_ due to scattering at _p_
            float pdf;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
                Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                        pdf;
            }
			
			// From photonmap.cpp (using our LVolumePhoton)
			ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);
			Lv += LVolumePhoton(volumeMap, nVolumePaths, nLookup, lookupBuf, vr, p, w, maxDistSquared);

        }
        ++sampOffset;
    }
    *T = Tr;
    return Lv * step * causticScale;
}

/* Copied from single.cpp */
Spectrum VolumePhotonIntegrator::Transmittance(const Scene *scene, 
		const Renderer *, const RayDifferential &ray, 
		const Sample *sample, RNG &rng, MemoryArena &arena) const
{
	if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

/* Modified from photonmap.cpp */
void VolumePhotonIntegrator::Preprocess(const Scene *scene, 
		const Camera *camera, const Renderer *renderer)
{
	if (scene->lights.size() == 0) return;
    // Declare shared variables for photon shooting
    Mutex *mutex = Mutex::Create();
    int nDirectPaths = 0;
    vector<Photon> volumetricPhotons;
    bool abortTasks = false;
    volumetricPhotons.reserve(nVolumePhotonsWanted);
    uint32_t nshot = 0;
    vector<Spectrum> rpReflectances, rpTransmittances;

    // Compute light power CDF for photon shooting
    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

    // Run parallel tasks for photon shooting
    ProgressReporter progress(nVolumePhotonsWanted, "Shooting volume photons");
    vector<Task *> photonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        photonShootingTasks.push_back(new VolumePhotonShootingTask(
            i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks, nDirectPaths,
            volumetricPhotons, rpReflectances, rpTransmittances, nshot, lightDistribution,
			scene, renderer, stepSize));
    EnqueueTasks(photonShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < photonShootingTasks.size(); ++i)
        delete photonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();

	if (volumetricPhotons.size() > 0)
        volumeMap = new KdTree<Photon>(volumetricPhotons);
}

#define RAY_EPSILON 0.00001

/* Modified from photonmap.cpp */
void VolumePhotonShootingTask::Run() {
    // Declare local variables for _PhotonShootingTask_
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localVolumePhotons;
    uint32_t totalPaths = 0;
    bool volumeDone = (integrator->nVolumePhotonsWanted == 0);
    PermutedHalton halton(6, rng);
    vector<Spectrum> localRpReflectances, localRpTransmittances;
	VolumeRegion *volumeRegion = scene->volumeRegion;
    while (true) {
        // Follow photon paths for a block of samples
        const uint32_t blockSize = 1024;
        for (uint32_t i = 0; i < blockSize; ++i) {
            float u[6];
            halton.Sample(++totalPaths, u);
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];

            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);
            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
            if (!alpha.IsBlack()) {
				
                // check for intersection with volume
                float t0, t1;
                if (!volumeRegion->IntersectP(photonRay, &t0, &t1)) continue;
                if (t0 > photonRay.mint + RAY_EPSILON)
                    photonRay = RayDifferential(photonRay(t0), photonRay.d, photonRay, RAY_EPSILON);
                
                bool scatter = false;
                int bounces = 0;
                
                // bounce photon
                while (volumeRegion->IntersectP(photonRay, &t0, &t1) &&
                       bounces < integrator->maxPhotonDepth + 1)
                {
                    photonRay.mint = t0;
                    photonRay.maxt = t0;
                    photonRay.d = Normalize(photonRay.d);
                    
                    // ray march
                    float rayt0 = t0, rayt1 = t1;
                    float sigmaTE = 0.f, dist = 0.f;
                    int count = 1;
                    float uRand = (float)rand()/(float)RAND_MAX;
                    while(rayt0 < rayt1){
                        sigmaTE += volumeRegion->sigma_t(photonRay(rayt0), -photonRay.d, stepSize).y();
                        dist = -1 * log(uRand) / (sigmaTE / float(count));
						if(rayt0 < dist)
                            break;
                        else {
                            rayt0 += stepSize;
                            count++;
                        }
                    }
                    dist *= stepSize;
                    
                    // break if point is previous intersection or outside volume
                    if (dist < RAY_EPSILON || dist > (t1 - t0))
                        break;
                    
                    photonRay.maxt += dist;
                    
                    // check for a surface intersection
                    Intersection isect;
                    if (scene->Intersect(photonRay, &isect)) {
                        // if intersection is close enough to current ray segment
                        if ((photonRay.o - isect.dg.p).Length() < t0 + dist) {
                            // Handle photon/surface intersection
                            BSDF *photonBSDF = isect.GetBSDF(photonRay, arena);
                            Vector wo = -photonRay.d;
                            Vector wi;
                            float pdf;
                            BxDFType flags;
                            Spectrum fr = photonBSDF->Sample_f(wo, &wi, BSDFSample(rng), &pdf, BSDF_ALL, &flags);
                            if (fr.IsBlack() || pdf == 0.f) break;
                            Spectrum anew = alpha * fr * AbsDot(wi, photonBSDF->dgShading.nn) / pdf;
                            
                            // Possibly terminate photon path with Russian roulette
                            float continueProb = min(1.f, anew.y() / alpha.y());
                            if (rng.RandomFloat() > continueProb)
                                break;
                            alpha = anew / continueProb;
                            
                            // update ray
                            photonRay = RayDifferential(isect.dg.p, wi, photonRay, isect.rayEpsilon);
                            bounces++;
                            
                            scatter = true;
                            continue;
                        }
                    }
                    
                    Point p = photonRay(t0 + dist);
                    if (scatter) {
                        Photon vp(p, alpha, -photonRay.d);
                        localVolumePhotons.push_back(vp);
                    }
                    
                    // get scattering properties
                    float sigmaT = volumeRegion->sigma_t(photonRay(t0), photonRay.d, stepSize).y();
					float sigmaS = volumeRegion->sigma_s(photonRay(t0), photonRay.d, stepSize).y();
					float albedo = sigmaS / sigmaT;
                    // Possibly terminate photon path with Russian roulette
					if (uRand < albedo) {
						float uRand1 = ((float)rand()) / (float)RAND_MAX;
						float uRand2 = ((float)rand()) / (float)RAND_MAX;
                        
                        // scatter photon
						Vector d = Normalize(SampleHG(photonRay.d, 0.f, uRand1, uRand2));
						photonRay = RayDifferential(p, d, photonRay, RAY_EPSILON);
					}
					else break;
                    
					bounces++;
					scatter = true;
                }
                
            }
            arena.FreeAll();
        }

        // Merge local photon data with data in _PhotonIntegrator_
        MutexLock lock(mutex);

        // Give up if we're not storing enough photons
        if (abortTasks)
            return;
        if (nshot > 500000 &&
            (unsuccessful(integrator->nVolumePhotonsWanted,
                          volumePhotons.size(), blockSize))) {
            Error("Unable to store enough volume photons.  Giving up.\n");
            volumePhotons.erase(volumePhotons.begin(), volumePhotons.end());
            abortTasks = true;
            return;
        }
        progress.Update(localVolumePhotons.size());
        nshot += blockSize;

        if (!volumeDone) {
            integrator->nVolumePaths += blockSize;
            for (uint32_t i = 0; i < localVolumePhotons.size(); ++i)
                volumePhotons.push_back(localVolumePhotons[i]);
            if (volumePhotons.size() >= integrator->nVolumePhotonsWanted)
                volumeDone = true;
            localVolumePhotons.erase(localVolumePhotons.begin(),
                                     localVolumePhotons.end());
        }

        // Exit task if enough photons have been found
        if (volumeDone)
            break;
    }
}

/* Modified from photonmap.cpp */
VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params) {
    int nVolume = params.FindOneInt("volumephotons", 20000);
    int nUsed = params.FindOneInt("nused", 50);
    if (PbrtOptions.quickRender) nVolume = nVolume / 10;
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    bool finalGather = params.FindOneBool("finalgather", false);
    float maxDist = params.FindOneFloat("maxdist", .1f);
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
    float causticScale = params.FindOneFloat("causticscale", 1.f);
	float stepSize = params.FindOneFloat("stepsize", 0.1f);
    return new VolumePhotonIntegrator(nVolume, nUsed, maxSpecularDepth, maxPhotonDepth, maxDist,
                                      finalGather, gatherAngle, stepSize, causticScale);
}

