// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn
// 
// monofonIC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// monofonIC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
#pragma once

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cstdint>
#include <numeric>
#include <vector>
#include <general.hh>

namespace particle{
	
class container
{
public:
	std::vector<float > positions32_, velocities32_, mass32_;
	std::vector<double> positions64_, velocities64_, mass64_;
	
	std::vector<uint32_t> ids32_;
	std::vector<uint64_t> ids64_;
	
	bool bhas_individual_masses_;

	container() : bhas_individual_masses_(false) { }

	container(const container &) = delete;

	void allocate(size_t nump, bool b64reals, bool b64ids, bool bindividualmasses)
	{
		bhas_individual_masses_ = bindividualmasses;

		if( b64reals ){
			positions64_.resize(3 * nump);
			velocities64_.resize(3 * nump);
			positions32_.clear();
			velocities32_.clear();
			if( bindividualmasses ){
				mass64_.resize(nump);
				mass32_.clear();
			}
		}else{
			positions32_.resize(3 * nump);
			velocities32_.resize(3 * nump);
			positions64_.clear();
			velocities64_.clear();
			if( bindividualmasses ){
				mass32_.resize(nump);
				mass64_.clear();
			}
		}

		if( b64ids ){
			ids64_.resize(nump);
			ids32_.clear();
		}else{
			ids32_.resize(nump);
			ids64_.clear();
		}
	}

	const void* get_pos32_ptr() const{
		return reinterpret_cast<const void*>( &positions32_[0] );
	}

	void set_pos32(size_t ipart, size_t idim, float p){
		positions32_[3 * ipart + idim] = p;
	}

	const void* get_pos64_ptr() const{
		return reinterpret_cast<const void*>( &positions64_[0] );
	}

	inline void set_pos64(size_t ipart, size_t idim, double p){
		positions64_[3 * ipart + idim] = p;
	}

	inline const void* get_vel32_ptr() const{
		return reinterpret_cast<const void*>( &velocities32_[0] );
	}
	
	inline void set_vel32(size_t ipart, size_t idim, float p){
		velocities32_[3 * ipart + idim] = p;
	}

	const void* get_vel64_ptr() const{
		return reinterpret_cast<const void*>( &velocities64_[0] );
	}

	inline void set_vel64(size_t ipart, size_t idim, double p){
		velocities64_[3 * ipart + idim] = p;
	}

	const void* get_ids32_ptr() const{
		return reinterpret_cast<const void*>( &ids32_[0] );
	}

	void set_id32(size_t ipart, uint32_t id){
		ids32_[ipart] = id;
	}

	const void* get_ids64_ptr() const{
		return reinterpret_cast<const void*>( &ids64_[0] );
	}

	void set_id64(size_t ipart, uint64_t id){
		ids64_[ipart] = id;
	}

	const void* get_mass32_ptr() const{
		return reinterpret_cast<const void*>( &mass32_[0] );
	}

	void set_mass32(size_t ipart, float m){
		mass32_[ipart] = m;
	}

	const void* get_mass64_ptr() const{
		return reinterpret_cast<const void*>( &mass64_[0] );
	}

	void set_mass64(size_t ipart, double m){
		mass64_[ipart] = m;
	}

	size_t get_local_num_particles(void) const
	{
		return std::max(ids32_.size(),ids64_.size());
	}

	size_t get_global_num_particles(void) const
	{
		size_t local_nump = this->get_local_num_particles(), global_nump;
#ifdef USE_MPI
		MPI_Allreduce(reinterpret_cast<void *>(&local_nump), reinterpret_cast<void *>(&global_nump), 1,
					  MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
		global_nump = local_nump;
#endif
		return global_nump;
	}

	size_t get_local_offset( void ) const
	{
		size_t this_offset = 0;

		#ifdef USE_MPI
			int mpi_size, mpi_rank;
			MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
			MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
			
			std::vector<size_t> nump_p_task(mpi_size,0), off_p_task;
			size_t num_p_this_task = this->get_local_num_particles();
			MPI_Allgather( reinterpret_cast<const void*>(&num_p_this_task), 1, MPI_UNSIGNED_LONG_LONG,
				reinterpret_cast<void*>(&nump_p_task[0]), 1, MPI_UNSIGNED_LONG_LONG, MPI_COMM_WORLD );

			off_p_task.push_back( 0 );
			std::partial_sum(nump_p_task.begin(), nump_p_task.end(), std::back_inserter(off_p_task) );
			this_offset = nump_p_task.at(mpi_rank);
		#endif

		return this_offset;
	}

	void dump(void)
	{
		/*for (size_t i = 0; i < ids_.size(); ++i)
		{
			std::cout << positions_[3 * i + 0] << " " << positions_[3 * i + 1] << " " << positions_[3 * i + 2] << " "
					  << velocities_[3 * i + 0] << " " << velocities_[3 * i + 1] << " " << velocities_[3 * i + 2] << std::endl;
		}*/
	}
};

}