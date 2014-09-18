/***************************************************************************
                          AllocateFields.hh
                         -------------------
    begin                : Mon Oct 4 2010
    copyright            : (C) 2010 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef ALLOCATEFIELDS_HH
#define ALLOCATEFIELDS_HH

#include "../defs.hh"

class AllocateFields {
public:
    AllocateFields(const Idx_t & Nx,
                   const Idx_t & Ny);

    ~AllocateFields();

    void execute(VField_Edge_t* & EFD,
                 VField_Cell_t* & HFD,
                 VField_Edge_t* & JFD,
                 VField_t* & AlphaE,
                 VField_t* & AlphaH,
                 Inform & msg);

private:

    int _Nx;
    int _Ny;

    BConds<Vector_t, DIM, UniformCartesian<DIM>, Edge> _vbc_edge;

    VField_Edge_t* _EFD;
    VField_Cell_t* _HFD;
    VField_Edge_t* _JFD;
    VField_t* _alpha_e;
    VField_t* _alpha_h;

    FieldLayout_t *_FL;
    FieldLayout_Edge_t *_FL_edge;
    FieldLayout_Cell_t *_FL_cell;
    Mesh_t *_mesh;
    Mesh_t *_mesh_cc;
};

#endif
