//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PiecewiseLinearInterpolationMaterial.h"
#include "PiecewiseTabularBase.h"
//#include "Piecewise.h"
#include "DelimitedFileReader.h"
// MOOSE includes
#include "MooseVariableFE.h"

registerMooseObject("MooseApp", PiecewiseLinearInterpolationMaterial);

defineLegacyParams(PiecewiseLinearInterpolationMaterial);

InputParameters
PiecewiseLinearInterpolationMaterial::validParams()
{

  InputParameters params = Material::validParams();
  params.addClassDescription("Compute a property using a piecewise linear interpolation to define "
                             "its dependence on a variable");
  params.addRequiredParam<std::string>("property",
                                       "The name of the property this material will compute");
  params.addRequiredCoupledVar(
      "variable",
      "The name of the variable whose value is used as the abscissa in the interpolation");
  params.addParam<std::vector<Real>>("x", "The abscissa values");
  params.addParam<std::vector<Real>>("y", "The ordinate values");
  params.addParam<std::vector<Real>>("xy_data",
                                     "All function data, supplied in abscissa, ordinate pairs");
  params.addParam<Real>("scale_factor", 1.0, "Scale factor to be applied to the ordinate values");
  params.addParam<bool>(
      "extrapolation",
      false,
      "Use linear extrapolation to evaluate points that lie outside given data set domain. ");
  params.declareControllable("scale_factor");
  return params;
}

PiecewiseLinearInterpolationMaterial::PiecewiseLinearInterpolationMaterial(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _prop_name(getParam<std::string>("property")),
    _coupled_var(coupledValue("variable")),
    _scale_factor(getParam<Real>("scale_factor")),
    _extrap(getParam<bool>("extrapolation")),
    _property(declareProperty<Real>(_prop_name)),
    _dproperty(declarePropertyDerivative<Real>(_prop_name, getVar("variable", 0)->name()))
{
  std::vector<Real> x;
  std::vector<Real> y;
  std::pair<std::vector<Real>, std::vector<Real>> xy;

  if ((parameters.isParamValid("x")) || (parameters.isParamValid("y")))
  {
    if (!((parameters.isParamValid("x")) && (parameters.isParamValid("y"))))
      mooseError("In PiecewiseLinearInterpolationMaterial ",
                 _name,
                 ": Both 'x' and 'y' must be specified if either one is specified.");

    if (parameters.isParamValid("xy_data"))
      mooseError("In PiecewiseLinearInterpolationMaterial ",
                 _name,
                 ": Cannot specify 'x', 'y', and 'xy_data' together.");

    x = getParam<std::vector<Real>>("x");
    y = getParam<std::vector<Real>>("y");
  }
  else if (parameters.isParamValid("xy_data"))
  {
    std::vector<Real> xy = getParam<std::vector<Real>>("xy_data");
    unsigned int xy_size = xy.size();
    if (xy_size % 2 != 0)
      mooseError("In PiecewiseLinearInterpolationMaterial ",
                 _name,
                 ": Length of data provided in 'xy_data' must be a multiple of 2.");

    unsigned int x_size = xy_size / 2;
    x.reserve(x_size);
    y.reserve(x_size);
    for (unsigned int i = 0; i < xy_size / 2; ++i)
    {
      x.push_back(xy[i * 2]);
      y.push_back(xy[i * 2 + 1]);
    }
  }
  else if (parameters.isParamValid("data_file"))
  {
    xy = buildFromFile();   
    x = xy.first;
    y = xy.second;
  }

  try
  {
    _linear_interp = libmesh_make_unique<LinearInterpolation>(x, y, _extrap);
  }
  catch (std::domain_error & e)
  {
    mooseError("In PiecewiseLinearInterpolationMaterial ", _name, ": ", e.what());
  }
}

void
PiecewiseLinearInterpolationMaterial::computeQpProperties()
{
  _property[_qp] = _scale_factor * _linear_interp->sample(_coupled_var[_qp]);
  _dproperty[_qp] = _scale_factor * _linear_interp->sampleDerivative(_coupled_var[_qp]);
}


std::pair<std::vector<Real>, std::vector<Real>>
PiecewiseLinearInterpolationMaterial::buildFromFile()
{
  // Input parameters
  const FileName & data_file_name = getParam<FileName>("data_file");
  const MooseEnum & format = getParam<MooseEnum>("format");
  unsigned int x_index = getParam<unsigned int>("x_index_in_file");
  unsigned int y_index = getParam<unsigned int>("y_index_in_file");
  bool xy_only = getParam<bool>("xy_in_file_only");

  // Check that other forms of input are not set.
  if (isParamValid("x") || isParamValid("y") || isParamValid("xy_data"))
    mooseError("In Piecewise ",
               _name,
               ": Cannot specify 'data_file' and 'x', 'y', or 'xy_data' together.");

  if (x_index == y_index)
    mooseError("In Piecewise ",
               _name,
               ": 'x_index_in_file' and 'y_index_in_file' are set to the same value.");

  // Read the data from CSV file
  MooseUtils::DelimitedFileReader reader(data_file_name, &_communicator);
  reader.setFormatFlag(format.getEnum<MooseUtils::DelimitedFileReader::FormatFlag>());
  reader.setComment("#");
  reader.read();
  const std::vector<std::vector<double>> & data = reader.getData();

  // Check the data format
  if (x_index >= data.size())
    mooseError("In Piecewise ",
               _name,
               ": The 'x_index_in_file' is out-of-range of the available data in '",
               data_file_name,
               "', which contains ",
               data.size(),
               " ",
               format,
               " of data.");

  if (y_index >= data.size())
    mooseError("In Piecewise ",
               _name,
               ": The 'y_index_in_file' is out-of-range of the available data in '",
               data_file_name,
               "', which contains ",
               data.size(),
               " ",
               format,
               " of data.");

  if (data.size() > 2 && xy_only)
    mooseError("In Piecewise ",
               _name,
               ": Read more than two ",
               format,
               " of data from file '",
               data_file_name,
               "'.  Did you mean to use \"format = ",
               format == "columns" ? "rows" : "columns",
               "\" or set \"xy_in_file_only\" to false?");

  // Update the input vectors to contained the desired data
  return std::make_pair(reader.getData(x_index), reader.getData(y_index));
}

std::pair<std::vector<Real>, std::vector<Real>>
PiecewiseLinearInterpolationMaterial::buildFromXandY()
{
  if (!((isParamValid("x")) && (isParamValid("y"))))
    mooseError("In PiecewiseLinearInterpolationMaterial ",
                _name,
                ": Both 'x' and 'y' must be specified if either one is specified.");

  if (isParamValid("xy_data"))
    mooseError("In PiecewiseLinearInterpolationMaterial ",
                _name,
                ": Cannot specify 'x', 'y', and 'xy_data' together.");

  return std::make_pair(getParam<std::vector<Real>>("x"), getParam<std::vector<Real>>("y"));
}

std::pair<std::vector<Real>, std::vector<Real>>
PiecewiseLinearInterpolationMaterial::buildFromXY()
{
  std::vector<Real> xy = getParam<std::vector<Real>>("xy_data");
  unsigned int xy_size = xy.size();
  if (xy_size % 2 != 0)
    mooseError("In PiecewiseLinearInterpolationMaterial ",
                _name,
                ": Length of data provided in 'xy_data' must be a multiple of 2.");

  unsigned int x_size = xy_size / 2;
  std::vector<Real> x;
  std::vector<Real> y;
  // x.reserve(data_size);
  // y.reserve(data_size);
  x.reserve(x_size);
  y.reserve(x_size);  
  // for (unsigned int i = 0; i < xy_size; i += 2)
  // {
  //   x.push_back(xy[i]);
  //   y.push_back(xy[i + 1]);
  // }
  for (unsigned int i = 0; i < xy_size / 2; ++i)
  {
    x.push_back(xy[i * 2]);
    y.push_back(xy[i * 2 + 1]);
  }
  return std::make_pair(x, y);
}
