#if defined(_MSC_VER)
#pragma once
#endif

#ifndef VOLUMEPHOTON_H
#define VOLUMEPHOTON_H

#include "volume.h"
#include "integrator.h"
#include "kdtree.h"


struct Photon;
struct ClosePhoton;
struct PhotonProcess;

class VolumePhotonIntegrator : public VolumeIntegrator
{
public:
	VolumePhotonIntegrator(int nvol, int nLookup, int maxspecdepth, int maxphotondepth,
        float maxdist, bool finalGather, float ga, float stepSize, float causticScale);
	~VolumePhotonIntegrator();
    
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera,
                    const Renderer *renderer);
    
	Spectrum Transmittance(const Scene *scene, const Renderer *renderer,
            const RayDifferential &ray, const Sample *sample, RNG &rng,
            MemoryArena &arena) const;

private:
	friend class VolumePhotonShootingTask;	

	uint32_t nVolumePhotonsWanted, nLookup;
    float maxDistSquared;
	int maxSpecularDepth, maxPhotonDepth;
	bool finalGather;
	float cosGatherAngle;

	int nVolumePaths;
	KdTree<Photon> *volumeMap;
    
    float stepSize;
	int tauSampleOffset, scatterSampleOffset;
    // scales the contribution from caustic photons
	float causticScale;
};

VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params);


#endif //VOLUMEPHOTON_H
