
/*
 *
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Nicolas Bourbaki
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */
 

#ifndef H__OGF_TOPOPT_COMMANDS_MESH_GROB_TOP_OPT_COMMANDS__H
#define H__OGF_TOPOPT_COMMANDS_MESH_GROB_TOP_OPT_COMMANDS__H

#include <OGF/TopOpt/common/common.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>
#include <geogram/mesh/mesh.h>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <armadillo>
#include <unistd.h>
#include <sys/wait.h> 
//#include "/home/hussein/Programming/project/graphite3_1.6.5/GraphiteThree/geogram/src/lib/exploragram/optimal_transport/optimal_transport_2d.h"

typedef std::pair<unsigned int, unsigned int> segment;
typedef std::pair<double , double> point_2;


namespace OGF {

	typedef index_t Face;
	typedef index_t Vertex;

    gom_class TopOpt_API MeshGrobTopOptCommands : public MeshGrobCommands {
    public:
        MeshGrobTopOptCommands() ;
        virtual ~MeshGrobTopOptCommands() ;

    gom_slots:

        //   Doxygen comments are parsed and used by Gomgen to
        //     generate tooltips.
        //   In addition to standard Doxygen tags, the following 
	//     tags can be used:
	//
        //      \menu  indicate a menu relative to current menu
        //           (MeshGrobTopOptCommands), or an absolute menu (starting
        //           with a '/') to insert the command in existing
        //           menus (for instance /Surface/Remesh)
	//          
        //      \advanced  all subsequent parameters are in the
        //        advanced section of the command (displayed when
        //        clicking on it)
		void add_matter(OGF::MeshGrob *m);
		void extract_mesh(GEO::Mesh *m);
		void send_mesh(OGF::MeshGrob *m);	
		void get_RVD(OGF::MeshGrob *m, OGF::MeshGrob * points);
 
	private:
		std::set<index_t> submesh;
		std::map<index_t, index_t> face_index_to_index;	
		std::map<index_t, index_t> vertex_index_to_index;
		std::map<index_t, index_t> index_to_index_vertex;
		std::set<index_t> fixed_faces;
		std::set<segment> borders;
		const std::string dj_dq_path = "/home/hussein/Desktop/mfem-3.3.2/examples/dq.txt";
		const std::string command    = "-m /home/hussein/Programming/project/graphite3_1.6.5/bin/my_mesh.mesh"; 
		const double epsilon = 0.0001;
		const double delta = 0.1;
	//	GEO::Attribute<bool> matter;
		// for rebuilding the mesh
		std::map<point_2, std::set<index_t>> point_to_vertex;
//		std::map<point_2, index_t> point_to_index;
		GEO::Mesh * my_mesh;
	
		void print_mesh(std::ostringstream &o, const MeshFacets & facets, const MeshVertices & vertices); 
		bool is_matter_face(const GEO::MeshFacets & facets, const GEO::MeshVertices & vertices,  const index_t f);
		bool is_matter_vertex(const double * coord);
		
		bool is_fixed_vertex(const double * coord);
		bool is_fixed_face(const GEO::MeshFacets& facets, const GEO::MeshVertices & vertices,  const index_t f);
		bool is_fixed_border(const segment& s, const MeshVertices& vertices);
		bool is_neumann_border(const segment& s, const MeshVertices& vertices); 
		bool is_neumann_vertex(const double * coord); 

		void get_vertices_with_matter(const GEO::MeshFacets & facets);
		void get_fixed_faces(const MeshFacets & facets, const MeshVertices & vertices);
		void get_border(const GEO::MeshFacetCornersStore& facets_corners, const MeshFacets & facets);	


		// for computing the derivative of the boundary/centroid 
		
		std::vector<index_t> get_surrounded_facets(const MeshFacetCornersStore & facet_corners, const MeshFacets& facets, const Face f, const Vertex v);
		std::set<index_t> get_centroid(const std::vector<Face> &incident_faces, const MeshFacets& facets);
	    arma::mat get_A_inverse(const std::set<index_t> & centroid, const MeshVertices & points); 
		arma::mat get_del_B(); 
		arma::mat get_del_AQ(const double qx, const double qy);

		GEO::vector<double> recieve_dj_dq(const std::string file_name); 
	    arma::mat derive_q_x_internal(const double qx, const double qy, const arma::mat &dB, const std::set<index_t> & centroid, const MeshVertices & points); 
 		void compute_dj_dx(OGF::MeshGrob *m,  OGF::MeshGrob * points); 
		void execute_cmd(const char *cmd);
    };
}

#endif

