
/*
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
 *  Contact for this Plugin: Nicolas Bourbaki or Hussein Houdrouge - hah51@mail.aub.edu
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
 
/*
 * TODO: 
 *	- send exec signal to mfem i.e automate the execution of mfem. -- ep done 
 *	- Recieve the dj/dq -- done 
 *  - apply the derivatives -- done  
 *  - update_the mesh -- done
 *  - optimal transport part
 */


#include <OGF/TopOpt/commands/mesh_grob_top_opt_commands.h>
#include <OGF/scene_graph/commands/scene_graph_commands.h>
#include <OGF/mesh/commands/mesh_grob_shapes_commands.h>
#include <geogram/basic/memory.h>
#include <geogram/parameterization/mesh_segmentation.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/matrix.h>
#include <exploragram/optimal_transport/optimal_transport_2d.h>

namespace OGF {

    MeshGrobTopOptCommands::MeshGrobTopOptCommands() {}
    MeshGrobTopOptCommands::~MeshGrobTopOptCommands() {}
	
	// check if the vertex belong to the middle of the domain, if yes then it is a void;     
	bool MeshGrobTopOptCommands::is_matter_vertex(const double * coord) {
		return not ((coord[1] > 0.25 && coord[1] < 0.75) && ((coord[0] > 0.25 && coord[0] < 0.75) || (coord[0] > 1.25 && coord[0] < 1.75)));	
	}

	// if the face has more than some percentage a void vertex then it is a void face
	/*
 	* return true if the face is matter , and false otherwise 
 	* -- need to change the name of this function to is_matter_face;
 	*/
	bool MeshGrobTopOptCommands::is_matter_face(const GEO::MeshFacets& facets, const GEO::MeshVertices & vertices,  const index_t f) {
		index_t count = 0;
		index_t nb_of_vertices = facets.nb_vertices(f);
		for(index_t v = 0; v < nb_of_vertices; ++v) {
			index_t tr_v = facets.vertex(f, v);
			const double * coord = vertices.point_ptr(tr_v);
			if(is_matter_vertex(coord)) {
				count++;
			}
		}
		double ratio = (count / ((double ) nb_of_vertices));
		return ratio >= 0.5;	
	}
	
	// set the void and the matter faces
	void MeshGrobTopOptCommands::add_matter(OGF::MeshGrob * points) {				
		GEO::AttributesManager& vam = points->vertices.attributes();
		GEO::Attribute<bool> v_matter(vam, "m");	

		for(index_t v = 0; v < points->vertices.nb(); ++v) {
			const double * coord = points->vertices.point_ptr(v);
			v_matter[v] = is_matter_vertex(coord);
		}
	}


	void MeshGrobTopOptCommands::set_matter(OGF::MeshGrob *m, OGF::MeshGrob * points) {
		GEO::AttributesManager& am =  (m->facets).attributes();
		GEO::Attribute<bool> matter(am, "m");
		GEO::AttributesManager& vam = points->vertices.attributes();
		GEO::Attribute<bool> v_matter(vam, "m");
		GEO::AttributesManager& f_am =  (m->facets).attributes();
		GEO::Attribute<index_t> charts(f_am, "chart");
		for(index_t i = 0; i < m->facets.nb(); i++) matter[i] = v_matter[charts[i]]; 
		m->update();
	}
	
//***************************************************************************************************************///
	// store the matter faces is sunmesh fiels ; after adding the matter I extracted the mesh with matter 
	void MeshGrobTopOptCommands::extract_mesh(GEO::Mesh *m) {
//		submesh.clear();
		GEO::AttributesManager& am =  (m->facets).attributes();
		GEO::Attribute<bool> matter(am, "m");
		index_t index = 0;
//		std::cout << "the number of facets is: " << m->facets.nb() << '\n';
		for(index_t i = 0;  i < (m->facets).nb(); ++i) {
			if(matter[i]) { 
				face_index_to_index[i] = index;
				index++;
				submesh.insert(i);
			}
		}
	}

//***************************************************************************************************************//	
	/*
 	 * This functions checks if the vertex belongs to the Dirichelet Boundary / or elemenets 
 	 *
 	 */
	bool MeshGrobTopOptCommands::is_fixed_vertex(const double * coord) {
		return (coord[0] == 0.0  && ((coord[1] >= 0.6  && coord[1] <= 1) || (coord[1] <= 0.3 && coord[1] >= 0)));
	}

	/*
 	*  This function checks wither the face is fixed or not by having at least one fixed vertex 
 	*
 	*/
	bool MeshGrobTopOptCommands::is_fixed_face(const GEO::MeshFacets& facets, const GEO::MeshVertices & vertices,  const index_t f) {
		index_t nb_of_vertices = facets.nb_vertices(f);
		if(nb_of_vertices > 3) return false;
		for(index_t v = 0; v < nb_of_vertices; ++v) {
			index_t tr_v = facets.vertex(f, v);
			const double * coord = vertices.point_ptr(tr_v);							
			if(is_fixed_vertex(coord)) return true;;
		}
		return false;	
	}

	/**
 	* 	Gather all the fixed face in a set;
 	*/
	void MeshGrobTopOptCommands::get_fixed_faces(const MeshFacets & facets, const MeshVertices & vertices) {
		index_t i = 0;
		for(auto&& f: submesh) { 
			if(i == (submesh.size() - 2)) break;
			if(is_fixed_face(facets, vertices, f)) fixed_faces.insert(f);
			i++;
		}
	}

	/**
 	*  In this function, I want to store the vertices with matter in a map from the true indexing in the old mesh
 	*  to the indexing [1..nb_of_vertices with matter]
 	*
 	*/
	void MeshGrobTopOptCommands::get_vertices_with_matter(const GEO::MeshFacets & facets) {
		GEO::AttributesManager& am =  facets.attributes();
		GEO::Attribute<bool> matter(am, "m");
		index_t count = 0;
		for(index_t i = 0; i < facets.nb(); i++) {
			if(matter[i]) {
				index_t nb_vertices = facets.nb_vertices(i);
				for(index_t v = 0; v < nb_vertices; ++v) { 
					index_t tr = facets.vertex(i, v);
					if(vertex_index_to_index.find(tr) == vertex_index_to_index.end()) {
						vertex_index_to_index[tr] = count;
						index_to_index_vertex[count] = tr;
						count++;
					}
				}	
			}
		}
	}
	
//***************************************************************************************************************//	
	/*
 	* Insert the edges that separates matter and a void in a set; 
 	* 	
 	*/
	void MeshGrobTopOptCommands::get_border(const GEO::MeshFacetCornersStore & facets_corners, const MeshFacets & facets) {
		
		//for(index_t i = 0; submesh.size(); ++i) {
		//index_t f = submesh[i];
		for(auto&& f: submesh) {
			index_t nb_of_vertices = facets.nb_vertices(f);
			if(nb_of_vertices > 10) continue;
			for(index_t j = 0; j < nb_of_vertices; ++j) {
				index_t c = facets.corner(f, j);
				index_t ad_face = facets_corners.adjacent_facet(c);
				if(submesh.find(ad_face) != submesh.end()) continue;
				// if we have matter and void 
				// we want to construct a segmet between this corner and the next one 
				index_t succ = facets.next_corner_around_facet(f, c);	
				index_t v1 = facets_corners.vertex(c);
				index_t v2 = facets_corners.vertex(succ);
				segment s1 = {v1, v2};
				segment s2 = {v2, v1};
				if((borders.find(s1) != borders.end()) || (borders.find(s2) != borders.end())) continue;
				borders.insert(s1);
			}
		}
	}

	/**
 	* maybe fix pointer issue 
 	*
 	*/
	bool MeshGrobTopOptCommands::is_fixed_border(const segment& s, const MeshVertices& vertices) {
		index_t v1 = s.first;
		index_t v2 = s.second;
		const double * coord_1 = vertices.point_ptr(v1);
		const double * coord_2 = vertices.point_ptr(v2);
		return  is_fixed_vertex(coord_1) && is_fixed_vertex(coord_2);
			
	} 
//***************************************************************************************************************//
	/**
 	* The neuman boundry is the boundry where the force is applied
 	*
 	*/
	bool MeshGrobTopOptCommands::is_neumann_vertex(const double * coord) {
		return (coord[0] == 2 && (coord[1] >= 0.2 && coord[1] <= 0.8));
	}

	bool MeshGrobTopOptCommands::is_neumann_border(const segment& s, const MeshVertices& vertices) {
		auto v1 = s.first;
		auto v2 = s.second;
		auto coord_1 = vertices.point_ptr(v1);
		auto coord_2 = vertices.point_ptr(v2);
		return is_neumann_vertex(coord_1) && is_neumann_vertex(coord_2);							
	}

//**************************************************************************************************************//
	void MeshGrobTopOptCommands::print_mesh(std::ostringstream & o, const MeshFacets & facets, const MeshVertices& vertices) {
		o << "MFEM mesh v1.0" << '\n';
	   	o << "#" << '\n' 
          << "# MFEM Geometry Types (see mesh/geom.hpp):" << '\n'
   		  << "#"  				 << '\n' 
		  << "# POINT       = 0" << '\n' 
   		  << "# SEGMENT     = 1" << '\n' 
		  << "# TRIANGLE    = 2" << '\n'
		  <<" # SQUARE      = 3" << '\n'
		  <<" # TETRAHEDRON = 4" << '\n'
		  <<" # CUBE        = 5" << '\n'
		  <<" #"   				 << '\n'; 
		  
		o << "dimension" << '\n';
		o << '2' << '\n';	
		o << '\n';
		o << "elements" << '\n';
		o << submesh.size() << '\n';
		//for(index_t i = 0; i < submesh.size(); ++i) {
		for(auto&& f : submesh) { 
		//	index_t f = submesh[i];
			if(fixed_faces.find(f) == fixed_faces.end()) o <<  '2' << ' ';
			else o << '1' << ' '; 

			o << '2' << ' ';
			// print the indeces of the vertices using my mapp
			index_t nb_of_vertices = facets.nb_vertices(f);			
			if(nb_of_vertices > 10) continue;
			for(index_t v = 0; v < nb_of_vertices; ++v) {
				index_t tr_v = facets.vertex(f, v);
				o << vertex_index_to_index[tr_v] << ' ';
			}
			o << '\n';
		}
		
		o << '\n';
		o << "boundary" << '\n';
		o << borders.size() << '\n';
		for(auto&& s: borders) {
			if(is_fixed_border(s, vertices)) o << '1' << ' ';
			else if(is_neumann_border(s, vertices)) o << '2' << ' ';
			else {
			  o << '3' << ' ';
			}	
			o << '1' << ' ' << vertex_index_to_index[s.first] << ' ' << vertex_index_to_index[s.second] << '\n'; 
		}
		
		o << '\n';
		o << "vertices" << '\n';
		o << vertex_index_to_index.size() << '\n';		
		o << '2' << '\n';
		for(auto&& e: index_to_index_vertex) {
	//		std::cout << "e fitst " << e.first << '\n';
			auto tr = e.second;
			const double * coord = vertices.point_ptr(tr);
			o << coord[0] << ' ' << coord[1] << '\n';
		} 
	}

	void MeshGrobTopOptCommands::send_mesh(OGF::MeshGrob *m) {
		(m->facets).connect();
		m->update();
//		m->save("original.mesh");
//		add_matter(m);
//		my_mesh->save("the_great_mesh.mesh");
//		my_mesh->update();
		extract_mesh(m);
	
		get_vertices_with_matter(m->facets);
		get_fixed_faces(m->facets, m->vertices);
		get_border(m->facet_corners, m->facets);

		std::ofstream out("/home/hussein/Desktop/mfem-3.3.2/examples/my_mesh.mesh");	
		std::ostringstream o;
		print_mesh(o, m->facets, m->vertices);		
		out << o.str() << '\n';
	}
//---------------------------------------------------------------------------------------------------------------

	std::vector<index_t> MeshGrobTopOptCommands::get_surrounded_facets(const MeshFacetCornersStore & facet_corners,const MeshFacets& facets, const Face f, const Vertex v) {
		std::vector<index_t> fs;
		Vertex local = facets.find_vertex(f, v);
		fs.push_back(f);
		index_t c = facets.corner(f, local);
		Face current = facet_corners.adjacent_facet(c);
		while (current != f) {
			if(current > facets.nb()) break;
			//std::cout << "the current face is " << current << std::endl;
			fs.push_back(current);
			local = facets.find_vertex(current, v);
			c = facets.corner(current, local);
			current = facet_corners.adjacent_facet(c);		
		}
		return fs;
	}

	std::set<index_t> MeshGrobTopOptCommands::get_centroid(const std::vector<Face> &incident_faces, const MeshFacets& facets) {
		std::set<Vertex> centroid; 
		if(incident_faces.size() == 0) {
			std::cout << "the number of incident faces is equal to 0"  << std::endl;
		}
		if(incident_faces.size() < 3) return centroid;
		GEO::AttributesManager& am = facets.attributes();
		GEO::Attribute<index_t> charts(am, "chart");
		for(index_t i = 0; i < incident_faces.size(); ++i) {
			Face f = incident_faces[i];
			centroid.insert(charts[f]);
		}
		// check if all of them are matter; 
		return centroid;
	}

	arma::mat MeshGrobTopOptCommands::get_A_inverse(const std::set<index_t> & centroid, const MeshVertices & points) {
		
		arma::mat A(2, 2);
	
		std::cout << "centroid number " << centroid.size() << std::endl;
		if(centroid.size() == 0) return A;
	
		auto it = centroid.begin();
		const double * x1 = points.point_ptr(*it);
		std::advance(it, 1);
		const double * x2 = points.point_ptr(*it);
		std::advance(it, 1);
		const double * x3 = points.point_ptr(*it);
	
		A(0, 0) = x2[0] - x1[0];
		A(0, 1) = x2[1] - x1[1];
		A(1, 0) = x3[0] - x1[0];
		A(1, 1) = x3[1] - x1[1];

		return inv(A);
	}

	arma::mat MeshGrobTopOptCommands::get_del_B() {
		arma::mat B(2, 9);
		// Row 1
		//----------------
		B(0, 0) =   -1;
		B(0, 1) =   -1;
		B(0, 2) =    1;
		B(0, 3) =    1;
		B(0, 4) =    0;
		B(0, 5) =    0;
		B(0, 6) =  0.5;
		B(0, 7) = -0.5;
		B(0, 8) = 	 0;
		// Row 2
		//----------------
		B(1, 0) =   -1;
		B(1, 1) =   -1;
		B(1, 2) =    0;
		B(1, 3) =    0;
		B(1, 4) =    1;
		B(1, 5) =    1;
		B(1, 6) =  0.5;
		B(1, 7) =    0;
		B(1, 8) = -0.5;			
		return B;
	}

	arma::mat MeshGrobTopOptCommands::get_del_AQ(const double qx, const double qy) {
		arma::mat AQ(2, 9);
		// Row 1
		// ---------------
		AQ(0, 0) = -qx;
		AQ(0, 1) = -qy;
		AQ(0, 2) =  qx;
		AQ(0, 3) =  qy;
		AQ(0, 4) =   0;
		AQ(0, 5) =   0;
		AQ(0, 6) =   0;
		AQ(0, 7) =   0;
		AQ(0, 8) =   0;
		// Row 2
		// ---------------
		AQ(1, 0) = -qx;
		AQ(1, 1) = -qy;
		AQ(1, 2) =   0;
		AQ(1, 3) =   0;
		AQ(1, 4) =  qx;
		AQ(1, 5) =  qy;
		AQ(1, 6) =   0;
		AQ(1, 7) =   0;
		AQ(1, 8) =   0;
  		return AQ;
	}

	arma::mat MeshGrobTopOptCommands::derive_q_x_internal(const double qx, const double qy, const arma::mat &dB, const std::set<index_t> & centroid, const MeshVertices & points) {
		auto A_inv = get_A_inverse(centroid, points);
		return  A_inv * (dB - 0.5 * get_del_AQ(qx, qy));
	}

	GEO::vector<double> MeshGrobTopOptCommands::recieve_dj_dq(const std::string file_name) {
		GEO::vector<double> dj_dq;
		std::fstream in(file_name);
		double x = 0.0;
		while((in >> x)) {dj_dq.push_back(x);}
		in.close();
		return dj_dq;
	}
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// separate this function into two but this ll add a storage cost 	
	void MeshGrobTopOptCommands::compute_dj_dx(OGF::MeshGrob *m,  OGF::MeshGrob * points, OGF::Mesh * omega) {

		// call mfem 
		execute_cmd();
		std::cout << "excuting mfem is done" << std::endl;	
		// create the vector
		// store the derivative of the centroid.	
		std::vector<double> dxx(points->vertices.nb());
		std::vector<double> dxy(points->vertices.nb());

	
		GEO::AttributesManager& am =  (m->facets).attributes();
		GEO::Attribute<index_t> charts(am, "chart");

		std::cout << "recieving dj dq" << std::endl;	
		auto dj_dq = recieve_dj_dq(dj_dq_path);
		std::cout << "done from recieving dj dq" << std::endl;	

		std::vector<bool> visited(dj_dq.size(), false); // we sould devide the size by two
		const auto del_B = get_del_B();
		std::cout << "the number of facets are " << m->facets.nb() << '\n';
		for(index_t i = 0; i < m->facets.nb(); ++i) {
			index_t nb_of_vertices = (m->facets).nb_vertices(i);
			for(Vertex v = 0; v < nb_of_vertices; ++v) {
				Vertex tr = m->facets.vertex(i, v);
				auto fem_index = vertex_index_to_index[tr];
				if(visited[fem_index]) continue;
				visited[fem_index] = true;
				auto incident_faces = get_surrounded_facets(m->facet_corners, m->facets, i, tr);
				auto centroid = get_centroid(incident_faces, m->facets);
				const double * coord = m->vertices.point_ptr(tr);
				// insert if statement //	
					auto dqdx = derive_q_x_internal(coord[0], coord[1], del_B, centroid, points->vertices);
					auto dqx = dj_dq[2 * fem_index];
			   		auto dqy = dj_dq[2 * fem_index + 1];
		
					dqdx.row(0) *= dqx;
					dqdx.row(1) *= dqy;
	
					unsigned int ii = 0;
					for(auto&& e: centroid) {
						dxx[e] = dqdx(0, ii) + dqdx(1, ii); 
						dxy[e] = dqdx(0, ii + 1) + dqdx(1, ii + 1);
						i += 2;
 					}
			} 
		}
		std::cout << "updating the mesh" << '\n';
 		update_mesh(dxx, dxy, points);	
		std::cout << "done from updating the mesh" << '\n';
		
//		const double * pts  = points->vertices.point_ptr(0);
//		std::vector<double> centroid(points->vertices.nb() * 3);
//		std::vector<double> w(points->vertices.nb());
//		compute_Laguerre_centroids_2d(omega, points->vertices.nb() , pts, centroid.data(), nullptr, false, 0, nullptr, 0, 0.0, nullptr, w.data(), 1000);
	}	

	void MeshGrobTopOptCommands::update_mesh(const std::vector<double> &dxx, const std::vector<double> &dxy, OGF::MeshGrob * points) {
		for(index_t i = 0; i < dxx.size(); ++i) {
			double * coord = points->vertices.point_ptr(i);
			coord[0] = (dxx[i] > 0.0 ? std::min(coord[0] + dxx[i], 2.0) : std::max(coord[0] + dxx[i], 0.0));
			coord[1] = (dxy[i] > 0.0 ? std::min(coord[1] + dxy[i], 1.0) : std::max(coord[1] + dxy[i], 0.0));
		}		
	}

//--------------------------------------------------------------------------------------
// Communication section 

	void MeshGrobTopOptCommands::execute_cmd() {
//		system("/home/hussein/Desktop/mfem-3.3.2/examples/ex2  -m /home/hussein/Desktop/mfem-3.3.2/examples/my_mesh.mesh");
		pid_t childpid = fork();
		if(childpid == 0) execlp("/home/hussein/Desktop/mfem-3.3.2/examples/ex2", "/home/hussein/Desktop/mfem-3.3.2/examples/ex2", "-m" ,"/home/hussein/Desktop/mfem-3.3.2/examples/my_mesh.mesh", NULL);
		else if(childpid > 0) wait(NULL);
	}
/*
	void MeshGrobTopOptCommands::get_RVD(OGF::MeshGrob *m, OGF::MeshGrob * points) {

		auto vv = recieve_dj_dq("/home/hussein/Desktop/mfem-3.3.2/examples/dq.txt");
		for(auto e : vv) {
			std::cout << "e in vv " << e << '\n';
		} 
		GEO::AttributesManager& am =  (m->facets).attributes();
		GEO::Attribute<index_t> charts(am, "chart");
		std::cout << "the charts size is: " << charts.size() << '\n';
		std::cout << "the number of faces is: " << (m->facets).nb() << '\n';
		for(index_t i = 0; i < m->facets.nb(); ++i) { 
			std::cout << charts[i] << '\n';
			const double * coord = (points->vertices).point_ptr(charts[i]);
			std::cout << coord[0] << ' ' << coord[1] << '\n';
			std::cout << (m->facets).nb_vertices(i) << '\n';
			index_t nb_of_vertices = (m->facets).nb_vertices(i);
			std::cout << nb_of_vertices << '\n';
			for(Vertex v = 0; v < nb_of_vertices; ++v) {
				Vertex tr = m->facets.vertex(i, v);
				auto incident_faces = get_surrounded_facets(m->facet_corners, m->facets, i, tr);
				std::cout << "the number of the incident faces: " << incident_faces.size() << '\n';
				auto centroid = get_centroid(incident_faces, m->facets);
				std::cout << "the number of centroid : " << centroid.size() << '\n';			
			} 
		}
	}*/
}

// replacing attribute by a function

