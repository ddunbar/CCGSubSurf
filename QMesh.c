#include <stdlib.h>
#include <string.h>

#include "CCGSubSurf.h"
#include "QMesh.h"

QMesh *qmesh_createGrid(int N) {
	QMesh *qm = malloc(sizeof(*qm));
	int i,j,idx;

	qm->numVerts = N*N;
	qm->numEdges = 2*N*(N-1);
	qm->numQuads = (N-1)*(N-1);

	qm->verts = malloc(sizeof(*qm->verts)*qm->numVerts);
	qm->edges = malloc(sizeof(*qm->edges)*qm->numEdges);
	qm->quads = malloc(sizeof(*qm->quads)*qm->numQuads);

	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			Vert *v = &qm->verts[j*N + i];

			v->co[0] = ((double) i - (N-1)/2.)/(N-1);
			v->co[1] = ((double) j - (N-1)/2.)/(N-1);
			v->co[2] = 0;
		}
	}

	idx = 0;
	for (i=0; i<N; i++) {
		for (j=0; j<N-1; j++) {
			Edge *e = &qm->edges[idx++];
			e->v[0] = &qm->verts[j*N + i];
			e->v[1] = &qm->verts[(j+1)*N + i];
		}
	}
	for (i=0; i<N-1; i++) {
		for (j=0; j<N; j++) {
			Edge *e = &qm->edges[idx++];
			e->v[0] = &qm->verts[j*N + i];
			e->v[1] = &qm->verts[j*N + i + 1];
		}
	}

	for (i=0; i<N-1; i++) {
		for (j=0; j<N-1; j++) {
			Quad *q = &qm->quads[j*(N-1) + i];
			q->v[0] = &qm->verts[(j+0)*N + (i+0)];
			q->v[1] = &qm->verts[(j+0)*N + (i+1)];
			q->v[2] = &qm->verts[(j+1)*N + (i+1)];
			q->v[3] = &qm->verts[(j+1)*N + (i+0)];
		}
	}

	return qm;
}

void qmesh_free(QMesh *qm) {
	free(qm->quads);
	free(qm->edges);
	free(qm->verts);
	free(qm);
}

///

static void qmesh_ccgMeshIFC_vertDataZero(CCGMeshHDL m, void *tv) {
	float *t = tv;
	t[0] = t[1] = t[2] = 0;
}
static int qmesh_ccgMeshIFC_vertDataEqual(CCGMeshHDL m, void *a, void *b) {
	return !memcmp(a, b, 12);
}
static void qmesh_ccgMeshIFC_vertDataCopy(CCGMeshHDL m, void *t, void *a) {
	memcpy(t, a, 12);
}
static void qmesh_ccgMeshIFC_vertDataAdd(CCGMeshHDL m, void *tav, void *bv) {
	float *ta = tav, *b = bv;

	ta[0] += b[0];
	ta[1] += b[1];
	ta[2] += b[2];
}
static void qmesh_ccgMeshIFC_vertDataSub(CCGMeshHDL m, void *tav, void *bv) {
	float *ta = tav, *b = bv;
	ta[0] -= b[0];
	ta[1] -= b[1];
	ta[2] -= b[2];
}
static void qmesh_ccgMeshIFC_vertDataMulN(CCGMeshHDL m, void *tav, double n) {
	float *ta = tav;
	ta[0] *= n;
	ta[1] *= n;
	ta[2] *= n;
}
static void qmesh_ccgMeshIFC_vertDataAvg4(CCGMeshHDL m, void *tv, void *av, void *bv, void *cv, void *dv) {
	float *t = tv, *a = av, *b = bv, *c = cv, *d = dv;
	t[0] = (a[0]+b[0]+c[0]+d[0])*.25;
	t[1] = (a[1]+b[1]+c[1]+d[1])*.25;
	t[2] = (a[2]+b[2]+c[2]+d[2])*.25;
}

void qmesh_syncCCGSubSurf(QMesh *qm, CCGSubSurf *ss) {
	int i;

	ccgSubSurf_initFullSync(ss);
	for (i=0; i<qm->numVerts; i++) {
		Vert *v = &qm->verts[i];

		ccgSubSurf_syncVert(ss, v, v->co);
	}
	for (i=0; i<qm->numEdges; i++) {
		Edge *e = &qm->edges[i];

		ccgSubSurf_syncEdge(ss, e, e->v[0], e->v[1]);
	}
	for (i=0; i<qm->numQuads; i++) {
		Quad *f = &qm->quads[i];

		ccgSubSurf_syncFace(ss, f, 4, (CCGVertHDL*) f->v);
	}
	ccgSubSurf_processSync(ss);
}

CCGSubSurf *qmesh_getCCGSubSurf(QMesh *qm, int levels) {
	CCGMeshIFC ifc = {0};
	CCGSubSurf *ss;

	ifc.vertUserSize = ifc.edgeUserSize = ifc.faceUserSize = 4;
	ifc.vertDataSize = 12;
	ifc.vertDataZero = qmesh_ccgMeshIFC_vertDataZero;
	ifc.vertDataEqual = qmesh_ccgMeshIFC_vertDataEqual;
	ifc.vertDataCopy = qmesh_ccgMeshIFC_vertDataCopy;
	ifc.vertDataAdd = qmesh_ccgMeshIFC_vertDataAdd;
	ifc.vertDataSub = qmesh_ccgMeshIFC_vertDataSub;
	ifc.vertDataMulN = qmesh_ccgMeshIFC_vertDataMulN;
	ifc.vertDataAvg4 = qmesh_ccgMeshIFC_vertDataAvg4;

	ss = ccgSubSurf_new(&ifc, qm, levels, NULL, NULL);

	ccgSubSurf_setUseAgeCounts(ss, 1, 0, 0, 0);

	return ss;
}
