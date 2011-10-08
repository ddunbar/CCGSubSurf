#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CCGSubSurf.h"
#include "QMesh.h"

int main(int argc, char** argv) {
  // Create the mesh and the subdivision surface.
  int subdivisionLevels = 5;
  QMesh *qm = qmesh_createGrid(21);
  CCGSubSurf *ss = qmesh_getCCGSubSurf(qm, subdivisionLevels);

  // Sync the mesh.,
  qmesh_syncCCGSubSurf(qm, ss);

  // Run our performance test, which is to do an N step lift of the mesh
  // into a monkey saddle. We want to do this in a way that also tests the
  // syncronization code, so we vary the lift rate from the exterior converges
  // (and is stable) earlier than the interior.
  //
  // We do two passes, one with and without full synchronization enabled, and we
  // report the approximate area under the final saddle (as a simple correctness
  // check).

  for (unsigned pass = 0; pass != 2; ++pass) {
    const unsigned numLiftSteps = 10;
    for (unsigned i = 0; i != numLiftSteps; ++i) {
      const float liftHeight = 1.0;
      float liftPercent = (float) i / (numLiftSteps - 1);
      qmesh_setToMonkeySaddleLift(qm, liftPercent, liftHeight,
                                  ss, /*applyAsFullSync=*/ pass==0);
    }

    // Compute the sum height of all final vertices.
    double areaUnderSaddle = 0;
    unsigned gridSize = ccgSubSurf_getGridSize(ss);
    CCGFaceIterator *fi = ccgSubSurf_getFaceIterator(ss);
    for (; !ccgFaceIterator_isStopped(fi); ccgFaceIterator_next(fi)) {
      CCGFace *f = ccgFaceIterator_getCurrent(fi);
      unsigned numVerts = ccgSubSurf_getFaceNumVerts(ss, f);
      for (unsigned S = 0; S != numVerts; ++S) {
        float (*faceGridData)[3] = ccgSubSurf_getFaceGridDataArray(ss, f, S);

        for (unsigned y = 0; y != gridSize - 1; ++y) {
          for (unsigned x = 0; x != gridSize - 1; ++x) {
            // We are lazy, and estimate the area under the quadrilateral as the
            // mean height times the area of the planar quadrilateral.
            float *v00 = faceGridData[(y+0)*gridSize + (x+0)];
            float *v01 = faceGridData[(y+0)*gridSize + (x+1)];
            float *v11 = faceGridData[(y+1)*gridSize + (x+1)];
            float *v10 = faceGridData[(y+1)*gridSize + (x+0)];
            float avgHeight = (v00[2] + v01[2] + v11[2] + v10[2]) / 4.0;
            float planarQuadArea =
              .5 * ((v11[0] - v00[0]) * (v01[1] - v10[1]) -
                    (v01[0] - v10[0]) * (v11[1] - v00[1]));
            areaUnderSaddle += planarQuadArea * fabs(avgHeight);
          }
        }
      }
    }
    printf("pass %d - area under saddle: %.6f\n", pass, areaUnderSaddle);
  }

  // Free the meshes.
  ccgSubSurf_free(ss);
  qmesh_free(qm);

  return 0;
}
