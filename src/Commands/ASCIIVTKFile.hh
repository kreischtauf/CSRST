/***************************************************************************
                           ASCIIVTKFile.hh
                         -------------------
    begin                : Sun Dec 4 2011
    copyright            : (C) 2011 by Christof Kraus
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

#ifndef ASCIIVTKFILE_HH
#define ASCIIVTKFILE_HH

#include "../defs.hh"
#include "../utils.hh"

class AsciiVtkFile {
public:
    AsciiVtkFile();

    template<class FT, class C>
    void addGridExtent(const Field<FT,DIM,Mesh_t,C> & exampleField);

    template<class C>
    void addScalarField(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField, const std::string & fieldName);
    template<class C>
    void addVectorField(const Field<Vector_t,DIM,Mesh_t,C> & vectorField, const std::string & fieldName);

    void writeFile(const std::string & baseFileName);

private:
    void addHeader();
    void addHeaderFileHeader();
    void addCoordinates(const Vector_t & dx);
    std::string addXCoordinates(const Vector_t & dx) const;
    std::string addYCoordinates(const Vector_t & dx) const;
    std::string addZCoordinates(const Vector_t & dx) const;
    void addGridExtent(const FieldLayout<DIM> & layout, const Vector_t & dx);
    void addScalarFieldHeader(const std::string & fieldName);
    void addVectorFieldHeader(const std::string & fieldName);
    void setLocalDomains(const FieldLayout<DIM> & layout);
    template<class C>
    std::string collectScalarData(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField);
    template<class C>
    std::string collectVectorData(const Field<Vector_t,DIM,Mesh_t,C> & vectorField);
    std::string concatenateFieldNames() const;
    std::string concatenateVectorFieldNames() const;
    std::string concatenateScalarFieldNames() const;
    std::string concatenateFields() const;
    std::string concatenateFieldHeaders() const;
    std::string concatenatePieceExtents(const std::string & baseFileName) const;
    void writeHeaderFile(const std::string & baseFileName);

    const double _cutRelativeValue;
    bool _localDomainsSet;
    bool _gridExtentAdded;
    bool _oneFieldAdded;
    bool _fileWritten;

    std::string _header;
    std::string _headerFileHeader;
    std::string _globalGridExtent;
    std::string _localGridExtent;
    std::string _coordinates;

    std::map<std::string, std::string> _fields;
    std::map<std::string, std::string> _fieldHeaders;
    std::vector<NDIndex<DIM> > _localDomains;

    Vector_t _origin;
};


inline
AsciiVtkFile::AsciiVtkFile():
    _cutRelativeValue(1e-6),
    _localDomainsSet(false),
    _gridExtentAdded(false),
    _oneFieldAdded(false),
    _fileWritten(false),
    _header(""),
    _headerFileHeader(""),
    _globalGridExtent(""),
    _localGridExtent(""),
    _coordinates(""),
    _fields(),
    _fieldHeaders(),
    _localDomains()
{
    addHeader();
}

template<class FT, class C>
void AsciiVtkFile::addGridExtent(const Field<FT,DIM,Mesh_t,C> & exampleField)
{
    FieldLayout<DIM> & layout = exampleField.getLayout();
    Vector_t dx(exampleField.get_mesh().get_meshSpacing(0),
                exampleField.get_mesh().get_meshSpacing(1));
    addGridExtent(layout, dx);
}

inline
void AsciiVtkFile::addGridExtent(const FieldLayout<DIM> & layout,
                                 const Vector_t & dx)
{
    if (!_gridExtentAdded) {
        if (!_localDomainsSet) {
            setLocalDomains(layout);
        }
        NDIndex<DIM> gDom = layout.getDomain();
        NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
        std::stringstream vtkPart;

        vtkPart << "\""
                << gDom[0].first() << " " << gDom[0].last() << " "
                << gDom[1].first() << " " << gDom[1].last() << " "
                << 1 << " " << 1
                << "\"";
        _globalGridExtent = vtkPart.str();
        vtkPart.str(std::string(""));
        vtkPart << "\""
                << localDomain[0].first() << " " << localDomain[0].last() << " "
                << localDomain[1].first() << " " << localDomain[1].last() << " "
                << 1 << " " << 1
                << "\"";
        _localGridExtent = vtkPart.str();

        _gridExtentAdded = true;

        addCoordinates(dx);
    }
}

inline
void AsciiVtkFile::setLocalDomains(const FieldLayout<DIM> & layout)
{
    Utils::getLocalDomains(layout, _localDomains);
    Utils::addGostCellToLocalDomains(layout, _localDomains);

    _localDomainsSet = true;
}

inline
void AsciiVtkFile::addScalarFieldHeader(const std::string & fieldName)
{
    std::stringstream vtkPart;

    vtkPart << indent_l0
            << "<PDataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"1\" format=\"ascii\"/>";

    _fieldHeaders.insert(std::pair<std::string, std::string>(fieldName, vtkPart.str()));
}

inline
void AsciiVtkFile::addVectorFieldHeader(const std::string & fieldName)
{
    std::stringstream vtkPart;

    vtkPart << indent_l0
            << "<PDataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"" << DIM << "\" format=\"ascii\"/>";

    _fieldHeaders.insert(std::pair<std::string, std::string>(fieldName, vtkPart.str()));
}

inline
void AsciiVtkFile::addHeader()
{
    std::stringstream vtkPart;
    short EndianTest_s;
    unsigned char *EndianTest = reinterpret_cast<unsigned char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;
    std::string endianness = (EndianTest_s == 1? "\"LittleEndian\"":"\"BigEndian\"");

    vtkPart << "<?xml version=\"1.0\"?>" << std::endl;
    vtkPart << indent_l0
            << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=" << endianness << ">\n";

    _header = vtkPart.str();

    addHeaderFileHeader();
}

inline
void AsciiVtkFile::addHeaderFileHeader()
{
    std::stringstream vtkPart;
    short EndianTest_s;
    unsigned char *EndianTest = reinterpret_cast<unsigned char*>(&EndianTest_s);
    EndianTest[0] = 1;
    EndianTest[1] = 0;
    std::string endianness = (EndianTest_s == 1? "\"LittleEndian\"":"\"BigEndian\"");

    vtkPart << "<?xml version=\"1.0\"?>" << std::endl;
    vtkPart << indent_l0
            << "<VTKFile type=\"PRectilinearGrid\" version=\"0.1\" byte_order=" << endianness << ">\n";

    _headerFileHeader = vtkPart.str();
}

inline
std::string AsciiVtkFile::addXCoordinates(const Vector_t & dx) const
{
    NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
    std::stringstream coordinatesStream;

    coordinatesStream << "<DataArray type=\"Float32\" name=\"X_COORDINATES\" NumberOfComponents=\"1\" format=\"ascii\">\n"
                      << indent_l1;

    for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
        coordinatesStream << i * dx[0] + _origin[0] << " ";
    }

    coordinatesStream << "\n"
                      << "</DataArray>";

    return coordinatesStream.str();
}

inline
std::string AsciiVtkFile::addYCoordinates(const Vector_t & dx) const
{
    NDIndex<DIM> localDomain = _localDomains[Ippl::myNode()];
    std::stringstream coordinatesStream;

    coordinatesStream << "<DataArray type=\"Float32\" name=\"Y_COORDINATES\" NumberOfComponents=\"1\" format=\"ascii\">\n"
                      << indent_l1;

    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        coordinatesStream <<  j * dx[1] + _origin[1] << " ";
    }

    coordinatesStream << "\n"
                      << "</DataArray>";

    return coordinatesStream.str();
}

inline
std::string AsciiVtkFile::addZCoordinates(const Vector_t & dx) const
{
    std::stringstream coordinatesStream;

    coordinatesStream << "<DataArray type=\"Float32\" name=\"Z_COORDINATES\" NumberOfComponents=\"1\" format=\"ascii\">\n"
                      << indent_l1 << 0.0 << "\n"
                      << "</DataArray>";
    return coordinatesStream.str();
}

inline
std::string AsciiVtkFile::concatenateFieldNames() const
{
    std::stringstream fieldNames;

    fieldNames << "Scalars=" << concatenateScalarFieldNames()
               << " Vectors=" << concatenateVectorFieldNames();
    return fieldNames.str();
}

template<class C>
void AsciiVtkFile::addScalarField(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField, const std::string & fieldName)
{
    _origin = scalarField.get_mesh().get_origin();
    if (!_gridExtentAdded) {
        addGridExtent(scalarField);
    }
    std::stringstream vtkPart;
    std::stringstream fieldNameAndComponents;

    fieldNameAndComponents << fieldName << ",1";
    vtkPart << indent_l0 << "<DataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"1\" format=\"ascii\">\n"
            << indent_l1 << collectScalarData(scalarField) << "\n"
            << indent_l0 << "</DataArray>";

    _fields.insert(std::pair<std::string, std::string>(fieldNameAndComponents.str(), vtkPart.str()));
    addScalarFieldHeader(fieldName);
    _oneFieldAdded = true;
}

template<class C>
void AsciiVtkFile::addVectorField(const Field<Vector_t,DIM,Mesh_t,C> & vectorField, const std::string & fieldName)
{
    _origin = vectorField.get_mesh().get_origin();
    if (!_gridExtentAdded) {
        addGridExtent(vectorField);
    }
    std::stringstream vtkPart;
    std::stringstream fieldNameAndComponents;

    fieldNameAndComponents << fieldName << ",3";
    vtkPart << indent_l0 << "<DataArray type=\"Float32\" Name=\"" << fieldName << "\" NumberOfComponents=\"" << DIM << "\" format=\"ascii\">\n"
            << indent_l1 << collectVectorData(vectorField) << "\n"
            << indent_l0 << "</DataArray>";

    _fields.insert(std::pair<std::string, std::string>(fieldNameAndComponents.str(), vtkPart.str()));
    addVectorFieldHeader(fieldName);
    _oneFieldAdded = true;
}

template<class C>
std::string AsciiVtkFile::collectScalarData(const Field<Scalar_t,DIM,Mesh_t,C> & scalarField)
{
    std::stringstream vtkPart;
    NDIndex<DIM> elem;
    const NDIndex<DIM> & localDomain = _localDomains[Ippl::myNode()];
    Scalar_t maxValue = max(abs(scalarField));

    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        elem[1]=Index(j,j);
        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
            elem[0]=Index(i,i);
            Scalar_t tmp = scalarField.localElement(elem);

            if (std::abs(tmp) / maxValue < _cutRelativeValue) {
                vtkPart << 0.0 << " ";
            } else {
                vtkPart << tmp << " ";
            }
        }
    }

    return vtkPart.str();
}

template<class C>
std::string AsciiVtkFile::collectVectorData(const Field<Vector_t,DIM,Mesh_t,C> & vectorField)
{
    std::stringstream vtkPart;
    NDIndex<DIM> elem;
    const NDIndex<DIM> & localDomain = _localDomains[Ippl::myNode()];
    Vector_t maxValue = max(vectorField);
    Vector_t minValue = min(vectorField);
    maxValue = Vector_t(std::max(maxValue(0),-minValue(0)), std::max(maxValue(1),-minValue(1)));

    for (int j = localDomain[1].first(); j <= localDomain[1].last(); ++ j) {
        elem[1]=Index(j,j);
        for (int i = localDomain[0].first(); i <= localDomain[0].last(); ++ i) {
            elem[0]=Index(i,i);
            Vector_t tmp = vectorField.localElement(elem);
            for (int l = 0; l < DIM; ++ l) {
                if (std::abs(tmp(l)) / maxValue(l) < _cutRelativeValue) {
                    vtkPart << 0.0 << " ";
                } else {
                    vtkPart << tmp(l) << " ";
                }
            }
        }
    }

    return vtkPart.str();
}

#endif
