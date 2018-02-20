#pragma once



#include <stdexcept>
#include <iostream>


#include <igl/readMESH.h>
#include <igl/faces_first.h>
#include <igl/boundary_facets.h>
#include <vector>
#include <GaussIncludes.h>
#include <ForceSpring.h>
#include <FEMIncludes.h>
#include <PhysicalSystemParticles.h>
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>

using namespace std;
using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring



typedef PhysicalSystemFEM<double, NeohookeanTet> NeohookeanTets;

typedef World<double,
	std::tuple<PhysicalSystemParticleSingle<double> *, NeohookeanTets *>,
	std::tuple<ForceSpringFEMParticle<double> *>,
	std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplictLinear<double, AssemblerEigenSparseMatrix<double>,
	AssemblerEigenVector<double> > MyTimeStepper;

namespace GAUSSDLL {

	class MyClass {
	public:
		MyClass(double grabbing_spring_stiffness, double spring_root_mass, double density, double youngsModulus, double poissonsRatio, double timestep, const char* mesh_path) {
			igl::readMESH(mesh_path, V, T, F); //"/meshesTetgen/Beam/Beam.mesh"
			igl::readMESH("C:\Users\cmurad\Documents\GAUSS\data\meshesTetgen\Beam\final_beams\4_.mesh", V1, T1, F1);

			//Use IGL's faces_first method
			//Pre-sort V, so Unity's mesh only needs V.topRows(F.maxCoeff()+1) to be rendered
			//Eigen::VectorXi IM;
		//	igl::faces_first(V, F, IM);
			//T = T.unaryExpr(bind1st(mem_fun( static_cast<Eigen::VectorXi::Scalar& (Eigen::VectorXi::*)(Eigen::VectorXi::Index)>(&Eigen::VectorXi::operator())), &IM)).eval();
			
			igl::boundary_facets(T, F);
			igl::boundary_facets(T1, F1);
			//Setup Physics
			
			tets = new NeohookeanTets(V, T);
			for (auto element : tets->getImpl().getElements()) {
				element->setDensity(density);
				element->setParameters(youngsModulus, poissonsRatio);
			}
			

			// // Pinned particle to attach spring for dragging
			pinned_point = new PhysicalSystemParticleSingle<double>();
			pinned_point->getImpl().setMass(spring_root_mass);
			auto fem_attached_pos = PosFEM<double>(&tets->getQ()[0], 0, &tets->getImpl().getV());

			spring_stiffness = grabbing_spring_stiffness; //4000.0;
			double spring_rest_length = 0.1;
			forceSpring = new ForceSpringFEMParticle<double>(fem_attached_pos, // TODO compare getV to V. Get rid of double use of index
				PosParticle<double>(&pinned_point->getQ()),
				spring_rest_length, 0.0);
			
			world.addSystem(pinned_point);
			world.addForce(forceSpring);
			world.addSystem(tets);
			fixDisplacementMin(world, tets);
			world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

			reset_world(world);

			N = createEmbeddedMesh(V1, tets, world);

			stepper = new MyTimeStepper(timestep);
			stepper->step(world);

			std::cout << F << std::endl;
			
		}

		void reset_world(MyWorld &world) {
			auto q = mapStateEigen(world); // TODO is this necessary?
			q.setZero();
		}



		int increment(bool grabbed, int index, double x, double y, double z) {
			
			if (grabbed) {
				auto q_part = mapDOFEigen(pinned_point->getQ(), world);
				Eigen::MatrixXd part_pos = q_part;

				//mouse down
				auto pinned_q = mapDOFEigen(pinned_point->getQ(), world);
				pinned_q = Eigen::RowVector3d(x, y, z);

				auto fem_attached_pos = PosFEM<double>(&tets->getQ()[index], index, &tets->getImpl().getV());
				forceSpring->getImpl().setPosition0(fem_attached_pos);
				forceSpring->getImpl().setStiffness(spring_stiffness);
			}
			else {
				forceSpring->getImpl().setStiffness(0.0);
			}


			stepper->step(world);
			return m_i;
		}

		Eigen::MatrixXd createEmbeddedMesh(Eigen::MatrixXd EmbeddingV, NeohookeanTets *tets, MyWorld &world) {
			Eigen::MatrixXd Nc;
			Nc.resize(3 * EmbeddingV.rows(), 3 * tets->getImpl().getV().rows()); //Embedded x originalMesh

			Nc.setZero();
			std::cout << Nc.rows() << ", " << Nc.cols() << std::endl << std::endl;
			// std::cout<<EmbeddingV<<std::endl<<std::endl;

			Eigen::VectorXd reachedThisShapeFunction(3 * tets->getImpl().getV().rows());

			for (unsigned int i = 0; i<EmbeddingV.rows(); ++i)
			{
				// std::cout<<"Row i "<<i<<std::endl;
				// std::cout<<EmbeddingV.row(i)<<std::endl;
				reachedThisShapeFunction.setZero();

				for (auto element : tets->getImpl().getElements())
				{
					int i0 = (element->EnergyKineticNonLumped::m_qDofs)[0]->getLocalId();
					int i1 = (element->EnergyKineticNonLumped::m_qDofs)[1]->getLocalId();
					int i2 = (element->EnergyKineticNonLumped::m_qDofs)[2]->getLocalId();
					int i3 = (element->EnergyKineticNonLumped::m_qDofs)[3]->getLocalId();

					double ro_array[3];
					Eigen::RowVectorXd::Map(ro_array, EmbeddingV.row(i).cols()) = EmbeddingV.row(i);
					double phi0 = element->EnergyKineticNonLumped::template phi<0>(ro_array);
					double phi1 = element->EnergyKineticNonLumped::template phi<1>(ro_array);
					double phi2 = element->EnergyKineticNonLumped::template phi<2>(ro_array);
					double phi3 = element->EnergyKineticNonLumped::template phi<3>(ro_array);
					if ((phi0<-1e-6 || phi0>(1 + 1e-6)) || (phi1<-1e-6 || phi1>(1 + 1e-6)) || (phi2<-1e-6 || phi2>(1 + 1e-6)) || (phi3<-1e-6 || phi3>(1 + 1e-6)))
					{
						// std::cout<<"phi val: "<<phi0<<", "<<phi1<<", "<<phi2<<", "<<phi3<<std::endl;
						phi0 = 0;
						phi1 = 0;
						phi2 = 0;
						phi3 = 0;
					}


					// std::cout<<"phi val: "<<phi0<<", "<<phi1<<", "<<phi2<<", "<<phi3<<std::endl;
					// std::cout<<"local Ids: "<<i0<<", "<<i1<<", "<<i2<<", "<<i3<<std::endl;
					// std::cout<<"DOFs: "<<
					// mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[0], world.getState()).transpose()<<", "<<
					// mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[1], world.getState()).transpose()<<", "<<
					// mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[2], world.getState()).transpose()<<", "<<
					// mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[3], world.getState()).transpose()<<
					// std::endl;

					// std::cout<<reachedThisShapeFunction.transpose()<<std::endl;
					if (reachedThisShapeFunction(i0)<1e-6) {
						// std::cout<<"i0 "<<i0<<std::endl;
						Nc(3 * i, i0) += phi0;
						Nc(3 * i + 1, i0 + 1) += phi0;
						Nc(3 * i + 2, i0 + 2) += phi0;
						reachedThisShapeFunction(i0) = phi0;
					}

					if (reachedThisShapeFunction(i1)<1e-6) {
						// std::cout<<"i1 "<<i1<<std::endl;
						Nc(3 * i, i1) += phi1;
						Nc(3 * i + 1, i1 + 1) += phi1;
						Nc(3 * i + 2, i1 + 2) += phi1;
						reachedThisShapeFunction(i1) = phi1;
					}

					if (reachedThisShapeFunction(i2)<1e-6) {
						// std::cout<<"i2 "<<i2<<std::endl;
						Nc(3 * i, i2) += phi2;
						Nc(3 * i + 1, i2 + 1) += phi2;
						Nc(3 * i + 2, i2 + 2) += phi2;
						reachedThisShapeFunction(i2) = phi2;
					}

					if (reachedThisShapeFunction(i3)<1e-6) {
						// std::cout<<"i3 "<<i3<<std::endl;
						Nc(3 * i, i3) += phi3;
						Nc(3 * i + 1, i3 + 1) += phi3;
						Nc(3 * i + 2, i3 + 2) += phi3;
						reachedThisShapeFunction(i3) = phi3;
					}
				}
				// std::cout<<"---------------------"<<std::endl;
				// std::cout<<std::endl;    
			}

			// std::cout<<Nc<<std::endl;
			return Nc;
		}

		Eigen::MatrixXd getCurrentVertPositions(MyWorld &world, NeohookeanTets *tets) {
			// Eigen::Map<Eigen::MatrixXd> q(mapStateEigen<0>(world).data(), V.cols(), V.rows()); // Get displacements only
			auto q = mapDOFEigen(tets->getQ(), world);
			Eigen::VectorXd q1 = N*q;
			Eigen::Map<Eigen::MatrixXd> dV(q1.data(), V1.cols(), V1.rows()); // Get displacements only

			return V1 + dV.transpose();
		}

		void get_mesh_sizes(int *n_face_indices, int *n_vertices) {
			*n_face_indices = F1.rows() * F1.cols();
			*n_vertices = V1.rows() * V1.cols();
		}

		void get_updated_mesh(double *verts_flat, int *face_indices_flat) {
			// Arrays must be preallocated
			Eigen::Map<Eigen::MatrixXd>(verts_flat, V1.cols(), V1.rows()) = getCurrentVertPositions(world, tets).transpose(); // swapped rows and cols?
			//verts_flat = V.transpose().data();
			Eigen::Map<Eigen::MatrixXi>(face_indices_flat, F1.cols(), F1.rows()) = F1.transpose();
			//face_indices_flat = F.data();
		}

	private:
		int m_i;

		Eigen::MatrixXd V; // Verts
		Eigen::MatrixXi T; // Tet indices
		Eigen::MatrixXi F; // Face indices

		Eigen::MatrixXd V1; // Verts
		Eigen::MatrixXi T1; // Tet indices
		Eigen::MatrixXi F1; // Face indices

		Eigen::MatrixXd N;

		NeohookeanTets *tets;
		MyWorld world;
		MyTimeStepper* stepper;
		PhysicalSystemParticleSingle<double> *pinned_point;
		ForceSpringFEMParticle<double> *forceSpring;

		double spring_stiffness;
	};

	MyClass* staticObject;
	extern "C" _declspec(dllexport) void* CreateMyObjectInstance(double grabbing_spring_stiffness, double spring_root_mass, double density, double youngsModulus, double poissonsRatio, double timestep, const char* mesh_path) {
		staticObject = new MyClass(grabbing_spring_stiffness, spring_root_mass, density, youngsModulus, poissonsRatio, timestep, mesh_path); // constructor is called
		return (void*) staticObject;

	}
	extern "C" _declspec(dllexport) int increment(void* object, bool grabbed, int index, double x, double y, double z) {
		return ((MyClass*)object)->increment(grabbed, index, x, y, z);
	}

	extern "C" _declspec(dllexport) void get_mesh_sizes(void* object, int *n_face_indices, int *n_vertices) {
		return ((MyClass*)object)->get_mesh_sizes(n_face_indices, n_vertices);
	}

	extern "C" _declspec(dllexport) void get_updated_mesh(void* object, double *verts_flat, int *face_indices_flat) {
		return ((MyClass*)object)->get_updated_mesh(verts_flat, face_indices_flat);
	}

	extern "C" _declspec(dllexport) void gauss_find_sizes(int * Vlen, int * Flen, char * node, char * ele);
	extern "C" _declspec(dllexport) void gauss_read_obj(double * V_flat, int * F_flat, int * Vlen, int * Flen, char * node, char * ele);
	extern "C" _declspec(dllexport) void gauss_step_sim();
}