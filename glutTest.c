#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "CCGSubSurf.h"

// Simple platform autodetection.
#if defined(__APPLE__)
#define PLATFORM_MACOSX
#elif defined(__linux__)
#define PLATFORM_LINUX
#elif defined(_WIN32)
#define PLATFORM_WIN32
#else
#error "Unknown platform"
#endif

#ifdef PLATFORM_MACOSX

#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>

#else

#ifdef PLATFORM_WIN32
#define WIN32_LEAN_AND_MEAN
#define WIN32_STRICT
#include <windows.h>
#endif

#include <GL/GL.h>
#include <GL/GLU.h>
#include <GLUT/glut.h>

#endif

///

#include "QMesh.h"

int intMax(int a, int b) {
  return (a<b)?b:a;
}

void ccgSubSurf_draw(CCGSubSurf *ss, int drawVerts, int drawEdges,
                     int drawInteriorEdges, int drawFaces, int useLighting) {
  int x, y, S;
  int edgeSize = ccgSubSurf_getEdgeSize(ss);
  int gridSize = ccgSubSurf_getGridSize(ss);

  if (drawVerts) {
    CCGVertIterator *vi= ccgSubSurf_getVertIterator(ss);

    glPointSize(4.0);
    glBegin(GL_POINTS);
    for (; !ccgVertIterator_isStopped(vi); ccgVertIterator_next(vi)) {
      CCGVert *v= ccgVertIterator_getCurrent(vi);
      glColor3ub(128, intMax(0, 255-ccgSubSurf_getVertAge(ss, v)*4), 0);
      glVertex3fv(ccgSubSurf_getVertData(ss, v));
    }
    glEnd();

    ccgVertIterator_free(vi);
  }

  if (drawEdges) {
    CCGEdgeIterator *ei= ccgSubSurf_getEdgeIterator(ss);

    for (; !ccgEdgeIterator_isStopped(ei); ccgEdgeIterator_next(ei)) {
      CCGEdge *e = ccgEdgeIterator_getCurrent(ei);
      float (*edgeData)[3] = ccgSubSurf_getEdgeDataArray(ss, e);
      glColor3ub(192, intMax(0, 255-ccgSubSurf_getEdgeAge(ss, e)*4), 0);
      glBegin(GL_LINE_STRIP);
      for (x=0; x<edgeSize; x++) {
        glVertex3fv(edgeData[x]);
        //glVertex3fv(ccgSubSurf_getEdgeData(ss, e, x));
      }
      glEnd();
    }

    ccgEdgeIterator_free(ei);
  }

  if (drawInteriorEdges) {
    CCGFaceIterator *fi= ccgSubSurf_getFaceIterator(ss);

    for (; !ccgFaceIterator_isStopped(fi); ccgFaceIterator_next(fi)) {
      CCGFace *f= ccgFaceIterator_getCurrent(fi);
      int numVerts= ccgSubSurf_getFaceNumVerts(ss, f);
      int ageCol = intMax(0,255 - ccgSubSurf_getFaceAge(ss, f)*4);

      for (S=0; S<numVerts; S++) {
        float (*faceGridData)[3] = ccgSubSurf_getFaceGridDataArray(ss, f, S);

        glColor3ub(0,ageCol,192);
        glBegin(GL_LINE_STRIP);
        for (x=0; x<gridSize; x++) {
          glVertex3fv(faceGridData[x]);
          //glVertex3fv(ccgSubSurf_getFaceGridData(ss, f, S, x, 0));
        }
        glEnd();
        glColor3ub(0,ageCol,0);
        for (y=1; y<gridSize-1; y++) {
          glBegin(GL_LINE_STRIP);
          for (x=0; x<gridSize; x++) {
            glVertex3fv(faceGridData[y*gridSize + x]);
            //glVertex3fv(ccgSubSurf_getFaceGridData(ss, f, S, x, y));
          }
          glEnd();
        }
        for (x=1; x<gridSize-1; x++) {
          glBegin(GL_LINE_STRIP);
          for (y=0; y<gridSize; y++) {
            glVertex3fv(faceGridData[y*gridSize + x]);
            //glVertex3fv(ccgSubSurf_getFaceGridData(ss, f, S, x, y));
          }
          glEnd();
        }
      }
    }

    ccgFaceIterator_free(fi);
  }

  if (drawFaces) {
    CCGFaceIterator *fi= ccgSubSurf_getFaceIterator(ss);

    glPolygonOffset(1.0, 1.0);
    glEnable(GL_POLYGON_OFFSET_FILL);
    if (useLighting) {
      glShadeModel(GL_FLAT);
      glEnable(GL_LIGHTING);
    }

    glColor3ub(128,128,128);
    for (; !ccgFaceIterator_isStopped(fi); ccgFaceIterator_next(fi)) {
      CCGFace *f= ccgFaceIterator_getCurrent(fi);
      int numVerts= ccgSubSurf_getFaceNumVerts(ss, f);

      for (S=0; S<numVerts; S++) {
        float (*faceGridData)[3] = ccgSubSurf_getFaceGridDataArray(ss, f, S);

        for (y=0; y<gridSize-1; y++) {
          float *c = NULL;

          glBegin(GL_QUAD_STRIP);
          for (x=0; x<gridSize; x++) {
            float *a = faceGridData[(y+0)*gridSize + x];
            float *b = faceGridData[(y+1)*gridSize + x];

            if (c && useLighting) {
              float a_bX = b[0]-a[0], a_bY = b[1]-a[1], a_bZ = b[2]-a[2];
              float a_cX = c[0]-a[0], a_cY = c[1]-a[1], a_cZ = c[2]-a[2];

              glNormal3f(       a_bY*a_cZ - a_bZ*a_cY,
                                a_bZ*a_cX - a_bX*a_cZ,
                                a_bX*a_cY - a_bY*a_cX);
            }

            glVertex3fv(a);
            glVertex3fv(b);

            c = b;
          }
          glEnd();
        }
      }
    }

    if (useLighting) {
      glDisable(GL_LIGHTING);
      glShadeModel(GL_SMOOTH);
    }
    glDisable(GL_POLYGON_OFFSET_FILL);

    ccgFaceIterator_free(fi);
  }
}

/****/

#define eSS_DrawVerts                   (1<<0)
#define eSS_DrawEdges                   (1<<1)
#define eSS_DrawInteriorEdges   (1<<2)
#define eSS_DrawFaces                   (1<<3)
#define eSS_UseLighting                 (1<<4)

typedef enum {
  kMutateKindOff = 0,
  kMutateKindPerturb
} MutationKind;

static double checkTime(void) {
  return (double) clock()/CLOCKS_PER_SEC;
}

void glutDrawString(float x, float y, char *string) {
  char c;

  glRasterPos2f(x, y);
  while ((c = *string++)) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
}

static QMesh *qmesh = NULL;
static CCGSubSurf *ss = NULL;
static int gGridSize = 9, gSubdivisionLevels = 2;
static int gSubsurfDrawFlags = eSS_DrawVerts|eSS_DrawEdges|eSS_UseLighting;
static int gMutationStep = 0, gNumMutationSteps = 10;
static MutationKind gMutate = kMutateKindPerturb;
static float gMutationVec[3], gMutationFactor = .05;
static int gRotate = 1, gFullSync = 1, gShowStatistics = 0;
static float gRotateTime = 0;
static Vert *gMutationVert;

void initMeshes(int gridSize, int subdivisionLevels) {
  while (1) {
    int faceCount = (gridSize-1)*(gridSize-1)*(1<<(2*subdivisionLevels));

    if (subdivisionLevels>1 && faceCount>10000000) {
      subdivisionLevels--;
    } else {
      break;
    }
  }

  if (!qmesh || gridSize != gGridSize) {
    if (qmesh) qmesh_free(qmesh);
    qmesh = qmesh_createGrid(gridSize);
    gMutationVert = NULL;
    if (ss) ccgSubSurf_free(ss);
    ss = NULL;
  }

  if (!ss) {
    ss = qmesh_getCCGSubSurf(qmesh, subdivisionLevels);
  } else {
    ccgSubSurf_setSubdivisionLevels(ss, subdivisionLevels);
  }
  qmesh_syncCCGSubSurf(qmesh, ss);

  gSubdivisionLevels = subdivisionLevels;
  gGridSize = gridSize;
}

#define kFPSFramesToTrack       64
static double gFPSFrameTimes[kFPSFramesToTrack];
static int gFPSCurFrame = 0;

void updateAndDrawStats(int height) {
  char buffer[512];
  int y = 0, textHeight = 15;
  double fps, lastTime, curTime, frameCount;

  if (gFPSCurFrame>=kFPSFramesToTrack) {
    lastTime = gFPSFrameTimes[gFPSCurFrame%kFPSFramesToTrack];
    frameCount = kFPSFramesToTrack;
  } else {
    lastTime = gFPSFrameTimes[0];
    frameCount = gFPSCurFrame;
  }

  gFPSFrameTimes[gFPSCurFrame%kFPSFramesToTrack] = curTime = checkTime();
  gFPSCurFrame++;

  glColor3f(1,1,1);

  glutDrawString(10, 10 + textHeight*y++, "(right-click for options)");
  y++;

  if (gShowStatistics) {
    sprintf(buffer, "Limit Mesh - # Verts: %d, # Edges: %d, # Faces: %d",
            ccgSubSurf_getNumFinalVerts(ss), ccgSubSurf_getNumFinalEdges(ss),
            ccgSubSurf_getNumFinalFaces(ss));
    glutDrawString(10, 10 + textHeight*y++, buffer);

    sprintf(buffer, "Grid Mesh - # Verts: %d, # Edges: %d, # Faces: %d",
            qmesh->numVerts, qmesh->numEdges, qmesh->numQuads);
    glutDrawString(10, 10 + textHeight*y++, buffer);

    sprintf(buffer, "Mutation: %s, Rate: %d, Factor: %2.2f", gMutate?"Yes":"No",
            gNumMutationSteps, gMutationFactor);
    glutDrawString(10, 10 + textHeight*y++, buffer);

    sprintf(buffer, "Grid Size: %d, Subdivision Levels: %d", gGridSize,
            gSubdivisionLevels);
    glutDrawString(10, 10 + textHeight*y++, buffer);

    sprintf(buffer, "Rotate: %s, Full Sync: %s", gRotate?"Yes":"No",
            gFullSync?"Yes":"No");
    glutDrawString(10, 10 + textHeight*y++, buffer);
  }

  fps = (double) frameCount/(curTime-lastTime);
  sprintf(buffer, "FPS: %2.2f", fps);
  glutDrawString(10, 10 + textHeight*y++, buffer);
}

void draw(void) {
  int w = glutGet(GLUT_WINDOW_WIDTH);
  int h =       glutGet(GLUT_WINDOW_HEIGHT);
  float lightPosition[4] = {0.0, 0.0, -1.0, 0.0};

  if (gRotate) {
    gRotateTime = checkTime();
  }

  glViewport(0, 0, w, h);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(35, (double) w/h, 0.1, 100.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,1);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

  glTranslatef(0, 0, -3);

  glRotated(gRotateTime*90, 1, 1, sin(gRotateTime));

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_POINT_SMOOTH);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);

  ccgSubSurf_draw(ss, gSubsurfDrawFlags&eSS_DrawVerts,
                  gSubsurfDrawFlags&eSS_DrawEdges,
                  gSubsurfDrawFlags&eSS_DrawInteriorEdges,
                  gSubsurfDrawFlags&eSS_DrawFaces,
                  gSubsurfDrawFlags&eSS_UseLighting);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, w, 0, h, -1, 1);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glFlush();

  updateAndDrawStats(h);

  glutSwapBuffers();

  glutPostRedisplay();
}

void timer(int value) {
  glutTimerFunc(1, timer, 0);

  switch (gMutate) {
  case kMutateKindOff:
    break;
  case kMutateKindPerturb:
    if (!gMutationVert || gMutationStep==0) {
      gMutationStep = 0;
      gMutationVert = &qmesh->verts[rand()%qmesh->numVerts];
      gMutationVec[0] = (1. - 2.*rand()/RAND_MAX);
      gMutationVec[1] = (1. - 2.*rand()/RAND_MAX);
      gMutationVec[2] = (1. - 2.*rand()/RAND_MAX);
    }

    gMutationVert->co[0] += gMutationFactor*gMutationVec[0]/gNumMutationSteps;
    gMutationVert->co[1] += gMutationFactor*gMutationVec[1]/gNumMutationSteps;
    gMutationVert->co[2] += gMutationFactor*gMutationVec[2]/gNumMutationSteps;

    gMutationStep = (gMutationStep+1)%gNumMutationSteps;

    if (gFullSync) {
      qmesh_syncCCGSubSurf(qmesh, ss);
    } else {
      ccgSubSurf_initPartialSync(ss);
      ccgSubSurf_syncVert(ss, (CCGVert*) gMutationVert, gMutationVert->co);
      ccgSubSurf_processSync(ss);
    }
    break;
  }
}

void keypress(unsigned char key, int x, int y) {
  if (key=='q' || key=='Q' || key==27) {
    exit(0);
  }
}

void menu_setMutateKind(int value) {
  gMutate = (MutationKind) value;
}

void menu_setMutateRate(int rate) {
  gNumMutationSteps = rate;
}

void menu_setMutateFactor(int factor) {
  gMutationFactor = factor/100.;
}

void menu_setSubdivisionLevels(int subdivisionLevels) {
  initMeshes(gGridSize, subdivisionLevels);
}

void menu_setGridSize(int gridSize) {
  initMeshes(gridSize, gSubdivisionLevels);
}

void menu_setSubsurfDisplayFlag(int bit) {
  gSubsurfDrawFlags ^= bit;
}

void menu_setMainOpts(int value) {
  if (value==0) {
    gRotate = !gRotate;
  } else if (value==1) {
    gFullSync = !gFullSync;
  } else if (value==2) {
    gShowStatistics = !gShowStatistics;
  }
}

void menuChange(int state) {
  if (state==GLUT_MENU_NOT_IN_USE) {
    gFPSCurFrame = 0;
  }
}

int main(int argc, char** argv) {
  initMeshes(gGridSize, gSubdivisionLevels);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowPosition(100, 100);
  glutInitWindowSize(640, 480);
  glutCreateWindow("CCGSubSurf GLUT Test");

  glutKeyboardFunc(keypress);
  glutDisplayFunc(draw);
  glutMenuStateFunc(menuChange);
  timer(0);

  {
    int mutate, mutateKind, mutateRate, mutateFactor;
    int subLevels, gridSize, subSurfDisp;

    mutateKind = glutCreateMenu(menu_setMutateKind);
    glutAddMenuEntry("Off", kMutateKindOff);
    glutAddMenuEntry("Randomly Perturb", kMutateKindPerturb);

    mutateRate = glutCreateMenu(menu_setMutateRate);
    glutAddMenuEntry("1", 1);
    glutAddMenuEntry("10", 10);
    glutAddMenuEntry("100", 100);
    glutAddMenuEntry("500", 500);
    glutAddMenuEntry("1000", 1000);

    mutateFactor = glutCreateMenu(menu_setMutateFactor);
    glutAddMenuEntry("0.01", 1);
    glutAddMenuEntry("0.05", 5);
    glutAddMenuEntry("0.10", 10);
    glutAddMenuEntry("0.20", 20);
    glutAddMenuEntry("0.50", 50);
    glutAddMenuEntry("1.00", 100);

    mutate = glutCreateMenu(0);
    glutAddSubMenu("Mutation Kind", mutateKind);
    glutAddSubMenu("Set Rate", mutateRate);
    glutAddSubMenu("Set Factor", mutateFactor);

    subLevels = glutCreateMenu(menu_setSubdivisionLevels);
    glutAddMenuEntry("1", 1);
    glutAddMenuEntry("2", 2);
    glutAddMenuEntry("3", 3);
    glutAddMenuEntry("4", 4);
    glutAddMenuEntry("5", 5);
    glutAddMenuEntry("6", 6);
    glutAddMenuEntry("7", 7);

    gridSize = glutCreateMenu(menu_setGridSize);
    glutAddMenuEntry("2", 2);
    glutAddMenuEntry("3", 3);
    glutAddMenuEntry("5", 5);
    glutAddMenuEntry("9", 9);
    glutAddMenuEntry("15", 15);
    glutAddMenuEntry("21", 21);
    glutAddMenuEntry("64", 64);
    glutAddMenuEntry("128", 128);

    subSurfDisp = glutCreateMenu(menu_setSubsurfDisplayFlag);

    glutAddMenuEntry("Draw Verts", eSS_DrawVerts);
    glutAddMenuEntry("Draw Edges", eSS_DrawEdges);
    glutAddMenuEntry("Draw Interior Edges", eSS_DrawInteriorEdges);
    glutAddMenuEntry("Draw Faces", eSS_DrawFaces);
    glutAddMenuEntry("Use Lighting", eSS_UseLighting);

    glutCreateMenu(menu_setMainOpts);
    glutAddMenuEntry("Rotate", 0);
    glutAddMenuEntry("Full Sync", 1);
    glutAddMenuEntry("Show Stats", 2);
    glutAddSubMenu("Mutate", mutate);
    glutAddSubMenu("SubSurf Display", subSurfDisp);
    glutAddSubMenu("Subdivision Levels", subLevels);
    glutAddSubMenu("Grid Size", gridSize);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
  }

  glutMainLoop();

  return 0;
}
