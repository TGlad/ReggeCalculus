#include <stdlib.h>
#include <glut.h>
#include <stdio.h>
#include "image.h"
#include "basics.h"
#include "SpaceTime.h"


static int width = 512;
static int height = 384;
static Image *image = NULL;

void init(void)
{
#if 0
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_FLAT);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  image = new Image(width, height);
#endif
}

void display(void)
{
#if 0
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      image->setPixel(0, 0, ScreenColour(Vector3(1, 1, 1)));
    }
  }
  glutSwapBuffers();
#endif
}

void update(void)
{
#if 0
  glutPostRedisplay();
#endif
}
void reshape(int w, int h)
{
  width = w;
  height = h;
}

void mouse(int button, int state, int x, int y) 
{
#if 0
   switch (button) {
      case GLUT_LEFT_BUTTON:
         if (state == GLUT_DOWN)
            glutIdleFunc(update); // starts the updating
         break;
      case GLUT_MIDDLE_BUTTON:
         if (state == GLUT_DOWN)
            glutIdleFunc(NULL);
         break;
      default:
         break;
   }
#endif
}

int main(int argc, char** argv)
{
#if 0
   glutInit(&argc, argv);   
   glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA);
   glutInitWindowSize(width, height); 
   glutInitWindowPosition(100, 100);
   glutCreateWindow(argv[0]);
   init();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutMouseFunc(mouse);
   glutMainLoop();
#else
  SpaceTime spacetime;
  spacetime.init();
#endif
   return 0;
}