/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */
#ifndef MRST_GEOMETRY_H_INCLUDED
#define MRST_GEOMETRY_H_INCLUDED

void compute_face_geometry(int ndims, const double *coords, int nfaces,
                           const unsigned int *nodepos, const int *facenodes,
                           double *fnormals, double *fcentroids,
                           double *fareas);
void compute_cell_geometry(int ndims,
                           const double* coords,
                           const unsigned* nodepos,
                           const int* facenodes,
                           const int* neighbors,
                           const double *fnormals,
                           const double *fcentroids,
                           int ncells, const unsigned int *facepos, const int *cellfaces,
                           double *ccentroids, double *cvolumes);

#endif /* MRST_GEOMETRY_H_INCLUDED */
