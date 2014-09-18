/***************************************************************************
                          AllocateFields.cpp
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

#include "AllocateFields.hh"

extern std::ofstream dbg;

AllocateFields::AllocateFields(const Idx_t & Nx,
                               const Idx_t & Ny):
    _Nx(Nx),
    _Ny(Ny),
    _EFD(NULL),
    _HFD(NULL),
    _JFD(NULL),
    _alpha_e(NULL),
    _alpha_h(NULL),
    _FL(NULL),
    _FL_edge(NULL),
    _FL_cell(NULL),
    _mesh(NULL)
{ }

AllocateFields::~AllocateFields()
{
    delete _EFD;
    _EFD = 0;
    delete _HFD;
    _HFD = 0;
    delete _JFD;
    _JFD = 0;
    delete _FL;
    _FL = 0;
    delete _mesh;
    _mesh = 0;
    delete _mesh_cc;
    _mesh_cc = 0;
    delete _alpha_e;
    _alpha_e = 0;
    delete _alpha_h;
    _alpha_h = 0;
}

void AllocateFields::execute(VField_Edge_t* & EFD,
                             VField_Cell_t* & HFD,
                             VField_Edge_t* & JFD,
                             VField_t* & AlphaE,
                             VField_t* & AlphaH,
                             Inform & msg)
{
    e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
    double spacing[] = {1.0, 1.0};

    for (int i = 0; i < 2 * DIM; ++ i) {
        _vbc_edge[i] = new ZeroGuardsAndZeroFace<Vector_t, DIM, UniformCartesian<DIM>, Edge>(i);
    }

    _mesh = new Mesh_t(Index(_Nx),
                       Index(_Ny),
                       spacing,
                       Vector_t(0.0));

    _mesh_cc = new Mesh_t(Index(_Nx + 1),
                          Index(_Ny + 1),
                          spacing,
                          Vector_t(0.0));

    _FL = new FieldLayout_t(*_mesh, decomp);
    _FL_edge = new FieldLayout_Edge_t(*_mesh, decomp);
    _FL_cell = new FieldLayout_Cell_t(*_mesh_cc, decomp);

    _EFD = new VField_Edge_t(*_mesh,
                             *_FL_edge,
                             // _vbcpp_edge,
                             GuardCellSizes<DIM>(GUARDCELLSIZE));
    *_EFD = 0.0;

    _HFD = new VField_Cell_t(*_mesh_cc,
                             *_FL_cell,
                             // _vbcpp_cell,
                             GuardCellSizes<DIM>(GUARDCELLSIZE));
    *_HFD = 0.0;

    _JFD = new VField_Edge_t(*_mesh,
                             *_FL_edge,
                             _vbc_edge,
                             GuardCellSizes<DIM>(GUARDCELLSIZE));
    *_JFD = 0.0;

    _alpha_e = new VField_t(*_mesh,
                            *_FL,
                            GuardCellSizes<DIM>(GUARDCELLSIZE));

    *_alpha_e = 0.0;

    _alpha_h = new VField_t(*_mesh,
                            *_FL,
                            GuardCellSizes<DIM>(GUARDCELLSIZE));

    *_alpha_h = 0.0;

    EFD = _EFD;
    HFD = _HFD;
    JFD = _JFD;
    AlphaE = _alpha_e;
    AlphaH = _alpha_h;
}
