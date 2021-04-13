/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#pragma once

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/brine.hh>
#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class HenryProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Henry { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

// Use a structured yasp grid
template<class TypeTag>
struct Grid<TypeTag, TTag::Henry> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Henry> { using type = HenryProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Henry>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::Brine<Scalar, Components::SimpleH2O<Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Henry>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = HenrySpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions to set salinity conveniently
template<class TypeTag>
struct UseMoles<TypeTag, TTag::Henry> { static constexpr bool value = false; };

} // end namespace Properties

/*!
 *
 */
template <class TypeTag>
class HenryProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    // copy pressure index for convenience
    enum { pressureIdx = Indices::pressureIdx };
    enum { dim = GridView::dimensionworld };

    //! The test is defined using mass fractions
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions!");

public:
    HenryProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();
        inflow_ = getParam<Scalar>("Problem.Inflow");
    }


    Scalar temperature() const
    { return 273.15 + 20; } // in [K]


    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        // Tipovi rubnih uvjeta dolaze ovdje
         return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
         PrimaryVariables values;
         // Dirichletov rubni uvjet dolazi ovdje
         return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
           NumEqVector values(0.0);
         // Neumannov rubni uvjet dolazi ovdje
           return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars;
         // Inicijalni uvjet dolazi ovdje
        return priVars;
    }


private:
    static constexpr Scalar eps_ = 1e-6;
    Scalar inflow_;
};

} // end namespace Dumux

