/* read_ppm.c
 * Reads in a PPM file with binary values and returns a 2D array of struct
 * ppm_pixel that stores RGB values
 * Jasmine Lei
 * 25 February 2022
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_ppm.h"

struct ppm_pixel** read_ppm(const char* filename, int* width, int* height) {
  FILE* file;
  char word[100];
  char temp[50];
  struct ppm_pixel ** pixels = NULL;

  file = fopen(filename, "rb");
  if (file == NULL) {
    return NULL;
  }

  // check for correct format
  fgets(word, 100, file);
  sscanf(word, "%s ", temp);
  if (strcmp(temp, "P6") != 0) {
    printf("Error: not an binary ppm file.\n");
    return NULL;
  }
  fgets(word, 100, file);
  // skip over comments
  while (word[0] == '#') {
    fgets(word, 100, file);
  }
  sscanf(word, "%d %d", width, height);
  // assuming max val is less than 256, and discard
  fgets(word, 100, file);

  // allocate space for 2D array of arrays
  pixels = malloc(sizeof(struct ppm_pixel *) * *height);
  if (pixels == NULL) {
    return NULL;
  }
  for (int i = 0; i < *height; i++) {
    pixels[i] = malloc(sizeof(struct ppm_pixel) * *width);
    if (pixels[i] == NULL) {
      return NULL;
    }
  }

  // populate array values
  for (int i = 0; i < *height; i++) {
    for (int j = 0; j < *width; j++) {
      fread(pixels[i][j].colors, 1, 3, file);
    }
  }

  fclose(file);
  return pixels;
}

extern void write_ppm(const char* filename, struct ppm_pixel** pxs, int w, int h) {
  FILE* file = fopen(filename, "wb");
  char text[100];

  if (file == NULL) {
    printf("Could not write to file.  Exiting...\n");
    exit(0);
  }

  // write header info
  fwrite("P6\n", 1, 3, file);
  text[0] = '\0';
  sprintf(text, "%d", w);
  fwrite(text, 1, strlen(text), file);
  fwrite(" ", 1, 1, file);
  text[0] = '\0';
  sprintf(text, "%d", h);
  fwrite(text, 1, strlen(text), file);
  fwrite("\n", 1, 1, file);
  text[0] = '\0';
  strcat(text, "255");
  fwrite(text, 1, strlen(text), file);
  fwrite("\n", 1, 1, file);

  // write from array to file, 3 chars at a time
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      fwrite(pxs[i][j].colors, 1, 3, file);
    }
  }
  fclose(file);
}
