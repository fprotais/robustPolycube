#pragma once
#include <ultimaille/all.h>
#include <array>


namespace rb_data_structure {
	struct Chart {
		int dim;
		int id;
		double value;
		std::vector<int> tetFacets;
	};
	struct Sorted_charts {
		std::array<int, 3> start_of_dim;
		std::vector<Chart> charts;
		std::vector<int> cf2chart;
		inline int size() const { return charts.size(); }
		inline const Chart& operator [](int i) const { return charts[i]; }
		inline Chart& operator [](int i) { return charts[i]; }

		inline int ranked(int dim, int num) const { return start_of_dim[dim] + num; }
		inline int dim_rank(int id) const { return charts[id].id - start_of_dim[charts[id].dim]; }
		inline int nb_charts(int dim) const { return (dim == 2 ? charts.size() - start_of_dim[dim] : start_of_dim[dim + 1] - start_of_dim[dim]); }
	};

	struct Block_decomposition {
		Block_decomposition() 
			: bnd_chart_id(m)
			, issuing_chart_id(m)
			, cell_volume(m)
			, block_coord(m)
			, float_coord(m)
			, int_coord(m)
		{};
		UM::Hexahedra m;
		UM::CellFacetAttribute<int> bnd_chart_id;
		UM::CellFacetAttribute<int> issuing_chart_id;
		UM::CellAttribute<double> cell_volume;
		UM::PointAttribute<UM::vec3> block_coord;
		UM::PointAttribute<UM::vec3> float_coord;
		UM::PointAttribute<UM::vec3> int_coord;

	};
	struct Final_mesh {
		Final_mesh()
			: original_bloc(m)
			, bnd_chart(m)
			, polycube_coord(m)
		{}
		UM::Hexahedra m;
		UM::CellAttribute<int> original_bloc;
		UM::CellFacetAttribute<int> bnd_chart;
		UM::PointAttribute<UM::vec3> polycube_coord;
	};
}

