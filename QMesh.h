#ifndef CCGSUBSURF_QMESH_H
#define CCGSUBSURF_QMESH_H

struct _CCGSubSurf;

typedef struct _Vert {
	float co[3];
} Vert;

typedef struct _Edge {
	Vert *v[2];
} Edge;

typedef struct _Quad {
	Vert *v[4];
} Quad;

typedef struct _QMesh {
	Vert *verts;
	Edge *edges;
	Quad *quads;
	int numVerts, numEdges, numQuads;
} QMesh;

QMesh *qmesh_createGrid(int N);
void qmesh_free(QMesh *qm);

CCGSubSurf *qmesh_getCCGSubSurf(QMesh *qm, int levels);

void qmesh_syncCCGSubSurf(QMesh *qm, CCGSubSurf *ss);

void qmesh_setToMonkeySaddleLift(
  QMesh *qm, float liftPercent, float liftHeight,
  struct _CCGSubSurf *ss, int applyAsFullSync);

#endif // CCGSUBSURF_QMESH_H
