/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#ifndef MRST_GEOMETRY_H_INCLUDED
#define MRST_GEOMETRY_H_INCLUDED

void compute_face_geometry(int ndims, double *coords, int nfaces,
                           unsigned *nodepos, int *facenodes,
                           double *fnormals, double *fcentroids,
                           double *fareas);
void compute_cell_geometry(int ndims, double *coords,
                           unsigned* nodepos, int *facenodes, int *neighbours,
                           double *fnormals,
                           double *fcentroids, int ncells,
                           unsigned* facepos, int *cellfaces,
                           double *ccentroids, double *cvolumes);

#endif /* MRST_GEOMETRY_H_INCLUDED */
