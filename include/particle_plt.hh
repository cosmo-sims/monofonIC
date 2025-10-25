// This file is part of monofonIC (MUSIC2)
// A software package to generate ICs for cosmological simulations
// Copyright (C) 2020 by Oliver Hahn & Bruno Marcos (this file)
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

//
// NOTE: part of this code, notably the translation into the first Brillouin 
//       zone (FBZ), was adapted from Fortran code written by Bruno Marcos
//
#pragma once

#include <general.hh>
#include <unistd.h> // for unlink

#include <iostream>
#include <fstream>

#include <random>
#include <map>

#include <cassert>

#include <particle_generator.hh>
#include <grid_fft.hh>
#include <math/mat3.hh>

#include <gsl/gsl_sf_hyperg.h>
inline double Hypergeometric2F1( double a, double b, double c, double x )
{
  return gsl_sf_hyperg_2F1( a, b, c, x);
}

#define PRODUCTION

namespace particle{
//! implement Joyce, Marcos et al. PLT calculation

class lattice_gradient{
private:
    const real_t boxlen_, aini_;
    const size_t ngmapto_, ngrid_, ngrid32_;
    const real_t mapratio_;//, XmL_;
    Grid_FFT<real_t,false> D_xx_, D_xy_, D_xz_, D_yy_, D_yz_, D_zz_;
    Grid_FFT<real_t,false> grad_x_, grad_y_, grad_z_;
    std::vector<vec3_t<real_t>> vectk_;
    std::vector<vec3_t<int>> ico_, vecitk_;

    bool is_even( int i ){ return (i%2)==0; }

    bool is_in( int i, int j, int k, const mat3_t<int>& M ){
        vec3_t<int> v({i,j,k});
        auto vv = M * v;
        return is_even(vv.x)&&is_even(vv.y)&&is_even(vv.z);
    }

    void init_D( lattice lattice_type )
    {
        constexpr real_t pi     = M_PI;
        constexpr real_t twopi  = 2.0*M_PI;
        constexpr real_t fourpi = 4.0*M_PI;
        const     real_t sqrtpi = std::sqrt(M_PI);
        const     real_t pi32   = std::pow(M_PI,1.5);

        //! === vectors, reciprocals and normals for the SC lattice ===
        const int charge_fac_sc = 1;
        const mat3_t<real_t> mat_bravais_sc{
            real_t{1.0}, real_t{0.0}, real_t{0.0},
            real_t{0.0}, real_t{1.0}, real_t{0.0},
            real_t{0.0}, real_t{0.0}, real_t{1.0}, 
        };
        const mat3_t<real_t> mat_reciprocal_sc{
            twopi, real_t{0.0}, real_t{0.0},
            real_t{0.0}, twopi, real_t{0.0},
            real_t{0.0}, real_t{0.0}, twopi,
        };
        const mat3_t<int> mat_invrecip_sc{
            2, 0, 0,
            0, 2, 0,
            0, 0, 2,
        };
        const std::vector<vec3_t<real_t>> normals_sc{
            {pi,real_t{0.},real_t{0.}},{-pi,real_t{0.},real_t{0.}},
            {real_t{0.},pi,real_t{0.}},{real_t{0.},-pi,real_t{0.}},
            {real_t{0.},real_t{0.},pi},{real_t{0.},real_t{0.},-pi},
        };
        

        //! === vectors, reciprocals and normals for the BCC lattice ===
        const int charge_fac_bcc = 2;
        const mat3_t<real_t> mat_bravais_bcc{
            real_t{1.0}, real_t{0.0}, real_t{0.5},
            real_t{0.0}, real_t{1.0}, real_t{0.5},
            real_t{0.0}, real_t{0.0}, real_t{0.5}, 
        };
        const mat3_t<real_t> mat_reciprocal_bcc{
            twopi, real_t{0.0}, real_t{0.0},
            real_t{0.0}, twopi, real_t{0.0},
            -twopi, -twopi, fourpi,
        };
        const mat3_t<int> mat_invrecip_bcc{
            2, 0, 0,
            0, 2, 0,
            1, 1, 1,
        };
        const std::vector<vec3_t<real_t>> normals_bcc{
            {real_t{0.0},pi,pi},{real_t{0.0},-pi,pi},{real_t{0.0},pi,-pi},{real_t{0.0},-pi,-pi},
            {pi,real_t{0.0},pi},{-pi,real_t{0.0},pi},{pi,real_t{0.0},-pi},{-pi,real_t{0.0},-pi},
            {pi,pi,real_t{0.0}},{-pi,pi,real_t{0.0}},{pi,-pi,real_t{0.0}},{-pi,-pi,real_t{0.0}}
        };
        

        //! === vectors, reciprocals and normals for the FCC lattice ===
        const int charge_fac_fcc = 4;
        const mat3_t<real_t> mat_bravais_fcc{
            real_t{0.0}, real_t{0.5}, real_t{0.0},
            real_t{0.5}, real_t{0.0}, real_t{1.0},
            real_t{0.5}, real_t{0.5}, real_t{0.0}, 
        };
        const mat3_t<real_t> mat_reciprocal_fcc{
            -fourpi, fourpi, twopi,
            real_t{0.0}, real_t{0.0}, twopi,
            fourpi, real_t{0.0}, -twopi,
        };
        const mat3_t<int> mat_invrecip_fcc{
            0, 1, 1,
            1, 0, 1,
            0, 2, 0,
        };
        const std::vector<vec3_t<real_t>> normals_fcc{
            {twopi,real_t{0.0},real_t{0.0}},{-twopi,real_t{0.0},real_t{0.0}},
            {real_t{0.0},twopi,real_t{0.0}},{real_t{0.0},-twopi,real_t{0.0}},
            {real_t{0.0},real_t{0.0},twopi},{real_t{0.0},real_t{0.0},-twopi},
            {+pi,+pi,+pi},{+pi,+pi,-pi},
            {+pi,-pi,+pi},{+pi,-pi,-pi},
            {-pi,+pi,+pi},{-pi,+pi,-pi},
            {-pi,-pi,+pi},{-pi,-pi,-pi},
        };
        
        //! select the properties for the chosen lattice
        const int ilat = lattice_type; // 0 = sc, 1 = bcc, 2 = fcc

        const auto mat_bravais     = (ilat==2)? mat_bravais_fcc : (ilat==1)? mat_bravais_bcc : mat_bravais_sc;
        const auto mat_reciprocal  = (ilat==2)? mat_reciprocal_fcc : (ilat==1)? mat_reciprocal_bcc : mat_reciprocal_sc;
        const auto mat_invrecip    = (ilat==2)? mat_invrecip_fcc : (ilat==1)? mat_invrecip_bcc : mat_invrecip_sc;
        const auto normals         = (ilat==2)? normals_fcc : (ilat==1)? normals_bcc : normals_sc;
        const auto charge_fac      = (ilat==2)? charge_fac_fcc : (ilat==1)? charge_fac_bcc : charge_fac_sc;

        const ptrdiff_t nlattice = ngrid_;
        const real_t dx = 1.0/real_t(nlattice);

        const real_t eta = 4.0; // Ewald cutoff shall be 4 cells
        const real_t alpha = 1.0/std::sqrt(2)/eta;
        const real_t alpha2 = alpha*alpha;
        const real_t alpha3 = alpha2*alpha;
        
        const real_t charge = 1.0/std::pow(real_t(nlattice),3)/charge_fac;
        const real_t fft_norm   = 1.0/std::pow(real_t(nlattice),3.0);
        const real_t fft_norm12 = 1.0/std::pow(real_t(nlattice),1.5);

        //! just a Kronecker \delta_ij
        auto kronecker = []( int i, int j ) -> real_t { return (i==j)? 1.0 : 0.0; };

        //! Ewald summation: short-range Green's function
        auto add_greensftide_sr = [&]( mat3_t<real_t>& D, const vec3_t<real_t>& d ) -> void {
            auto r = d.norm();
            if( r< 1e-14 ) return; // return zero for r=0

            const real_t r2(r*r), r3(r2*r), r5(r3*r2);
            const real_t K1( -alpha3/pi32 * std::exp(-alpha2*r2)/r2 );
            const real_t K2( (std::erfc(alpha*r) + 2.0*alpha/sqrtpi*std::exp(-alpha2*r2)*r)/fourpi );
            
            for( int mu=0; mu<3; ++mu ){
                for( int nu=mu; nu<3; ++nu ){
                    real_t dd( d[mu]*d[nu] * K1 + (kronecker(mu,nu)/r3 - 3.0 * (d[mu]*d[nu])/r5) * K2 );
                    D(mu,nu) += dd;
                    D(nu,mu) += (mu!=nu)? dd : 0.0;
                }
            }
        };

        //! Ewald summation: long-range Green's function
        auto add_greensftide_lr = [&]( mat3_t<real_t>& D, const vec3_t<real_t>& k, const vec3_t<real_t>& r ) -> void {
            real_t kmod2 = k.norm_squared();
            real_t term = std::exp(-kmod2/(4*alpha2))*std::cos(k.dot(r)) / kmod2 * fft_norm;
            for( int mu=0; mu<3; ++mu ){
                for( int nu=mu; nu<3; ++nu ){
                    auto dd = k[mu] * k[nu] * term;
                    D(mu,nu) += dd;
                    D(nu,mu) += (mu!=nu)? dd : 0.0;
                }
            }
        };

        //! checks if 'vec' is in the FBZ with FBZ normal vectors given in 'normals'
        auto check_FBZ = []( const auto& normals, const auto& vec ) -> bool {
            for( const auto& n : normals ){ 
                if( n.dot( vec ) > 1.0001 * n.dot(n) ){
                    return false;
                }
            }
            return true;
        };
        
        constexpr ptrdiff_t lnumber = 3, knumber = 3;
        const int numb = 1; //!< search radius when shifting vectors into FBZ

        vectk_.assign(D_xx_.memsize(),vec3_t<real_t>());
        ico_.assign(D_xx_.memsize(),vec3_t<int>());
        vecitk_.assign(D_xx_.memsize(),vec3_t<int>());

        #pragma omp parallel 
        {
            //... temporary to hold values of the dynamical matrix 
            mat3_t<real_t> matD(real_t(0.0));

            #pragma omp for
            for( ptrdiff_t i=0; i<nlattice; ++i ){
                for( ptrdiff_t j=0; j<nlattice; ++j ){
                    for( ptrdiff_t k=0; k<nlattice; ++k ){
                        // compute lattice site vector from (i,j,k) multiplying Bravais base matrix, and wrap back to box
                        const vec3_t<real_t> x_ijk({dx*real_t(i),dx*real_t(j),dx*real_t(k)});
                        const vec3_t<real_t> ar = (mat_bravais * x_ijk).wrap_abs();

                        //... zero temporary matrix
                        matD.zero();        

                        // add real-space part of dynamical matrix, periodic copies
                        for( ptrdiff_t ix=-lnumber; ix<=lnumber; ix++ ){
                            for( ptrdiff_t iy=-lnumber; iy<=lnumber; iy++ ){
                                for( ptrdiff_t iz=-lnumber; iz<=lnumber; iz++ ){      
                                    const vec3_t<real_t> n_ijk({real_t(ix),real_t(iy),real_t(iz)});            
                                    const vec3_t<real_t> dr(ar - mat_bravais * n_ijk);
                                    add_greensftide_sr(matD, dr);
                                }
                            }
                        }

                        // add k-space part of dynamical matrix
                        for( ptrdiff_t ix=-knumber; ix<=knumber; ix++ ){
                            for( ptrdiff_t iy=-knumber; iy<=knumber; iy++ ){
                                for( ptrdiff_t iz=-knumber; iz<=knumber; iz++ ){                      
                                    if(std::abs(ix)+std::abs(iy)+std::abs(iz) != 0){
                                        const vec3_t<real_t> k_ijk({real_t(ix)/nlattice,real_t(iy)/nlattice,real_t(iz)/nlattice});
                                        const vec3_t<real_t> ak( mat_reciprocal * k_ijk);

                                        add_greensftide_lr(matD, ak, ar );
                                    }
                                }
                            }
                        } 

                        D_xx_.relem(i,j,k) = matD(0,0) * charge;
                        D_xy_.relem(i,j,k) = matD(0,1) * charge;
                        D_xz_.relem(i,j,k) = matD(0,2) * charge;
                        D_yy_.relem(i,j,k) = matD(1,1) * charge;
                        D_yz_.relem(i,j,k) = matD(1,2) * charge;
                        D_zz_.relem(i,j,k) = matD(2,2) * charge;
                    }
                }
            }
        } // end omp parallel region

        // fix r=0 with background density (added later in Fourier space)
        D_xx_.relem(0,0,0) = 1.0/3.0;
        D_xy_.relem(0,0,0) = 0.0;
        D_xz_.relem(0,0,0) = 0.0;
        D_yy_.relem(0,0,0) = 1.0/3.0;
        D_yz_.relem(0,0,0) = 0.0;
        D_zz_.relem(0,0,0) = 1.0/3.0;

        D_xx_.FourierTransformForward();
        D_xy_.FourierTransformForward();
        D_xz_.FourierTransformForward();
        D_yy_.FourierTransformForward();
        D_yz_.FourierTransformForward();
        D_zz_.FourierTransformForward();

#ifndef PRODUCTION
        if (CONFIG::MPI_task_rank == 0)
            unlink("debug.hdf5");
        D_xx_.Write_to_HDF5("debug.hdf5","Dxx");
        D_xy_.Write_to_HDF5("debug.hdf5","Dxy");
        D_xz_.Write_to_HDF5("debug.hdf5","Dxz");
        D_yy_.Write_to_HDF5("debug.hdf5","Dyy");
        D_yz_.Write_to_HDF5("debug.hdf5","Dyz");
        D_zz_.Write_to_HDF5("debug.hdf5","Dzz");

        std::ofstream ofs2("test_brillouin.txt");
#endif
        using map_t = std::map<vec3_t<int>,size_t>;
        map_t iimap;
            
        //!=== Make temporary copies before resorting to std. Fourier grid ========!//
        Grid_FFT<real_t,false> 
            temp1({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
            temp2({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
            temp3({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0});

        temp1.FourierTransformForward(false);
        temp2.FourierTransformForward(false);
        temp3.FourierTransformForward(false);
            
        #pragma omp parallel for
        for( size_t i=0; i<D_xx_.size(0); i++ )
        {
            for( size_t j=0; j<D_xx_.size(1); j++ )
            {
                for( size_t k=0; k<D_xx_.size(2); k++ )
                {
                    temp1.kelem(i,j,k) = ccomplex_t(std::real(D_xx_.kelem(i,j,k)),std::real(D_xy_.kelem(i,j,k)));
                    temp2.kelem(i,j,k) = ccomplex_t(std::real(D_xz_.kelem(i,j,k)),std::real(D_yy_.kelem(i,j,k)));
                    temp3.kelem(i,j,k) = ccomplex_t(std::real(D_yz_.kelem(i,j,k)),std::real(D_zz_.kelem(i,j,k)));
                }
            }
        }
        D_xx_.zero(); D_xy_.zero(); D_xz_.zero();
        D_yy_.zero(); D_yz_.zero(); D_zz_.zero();

        
        //!=== Diagonalise and resort to std. Fourier grid ========!//
        #pragma omp parallel 
        {
            // thread private matrix representation
            mat3_t<real_t> D;
            vec3_t<real_t> eval, evec1, evec2, evec3_t;

            #pragma omp for
            for( size_t i=0; i<D_xx_.size(0); i++ )
            {
                for( size_t j=0; j<D_xx_.size(1); j++ )
                {
                    for( size_t k=0; k<D_xx_.size(2); k++ )
                    {
                        vec3_t<real_t> kv = D_xx_.get_k<real_t>(i,j,k);
                        
                        // put matrix elements into actual matrix
                        D(0,0) = std::real(temp1.kelem(i,j,k)) / fft_norm12;
                        D(0,1) = D(1,0) = std::imag(temp1.kelem(i,j,k)) / fft_norm12;
                        D(0,2) = D(2,0) = std::real(temp2.kelem(i,j,k)) / fft_norm12;
                        D(1,1) = std::imag(temp2.kelem(i,j,k)) / fft_norm12;
                        D(1,2) = D(2,1) = std::real(temp3.kelem(i,j,k)) / fft_norm12;
                        D(2,2) = std::imag(temp3.kelem(i,j,k)) / fft_norm12;

                        // compute eigenstructure of matrix
                        D.eigen(eval, evec1, evec2, evec3_t);
                        evec3_t /= (twopi*ngrid_);

                        // now determine to which modes on the regular lattice this contributes
                        vec3_t<real_t> ar = kv / (twopi*ngrid_);
                        vec3_t<real_t> a(mat_reciprocal * ar);
                        
                        // translate the k-vectors into the "candidate" FBZ
                        for( int l1=-numb; l1<=numb; ++l1 ){
                            for( int l2=-numb; l2<=numb; ++l2 ){
                                for( int l3=-numb; l3<=numb; ++l3 ){
                                    // need both halfs of Fourier space since we use real transforms
                                    for( int isign=0; isign<=1; ++isign ){
                                        const real_t sign = 2.0*real_t(isign)-1.0; 
                                        const vec3_t<real_t> vshift({real_t(l1),real_t(l2),real_t(l3)});

                                        vec3_t<real_t> vectk = sign * a + mat_reciprocal * vshift;

                                        if( check_FBZ( normals, vectk ) )
                                        {
                                            int ix = std::round(vectk.x*(ngrid_)/twopi);
                                            int iy = std::round(vectk.y*(ngrid_)/twopi);
                                            int iz = std::round(vectk.z*(ngrid_)/twopi);

                                            #pragma omp critical
                                            {iimap.insert( std::pair<vec3_t<int>,size_t>({ix,iy,iz}, D_xx_.get_idx(i,j,k)) );}

                                            temp1.kelem(i,j,k) = ccomplex_t(eval[2],eval[1]);
                                            temp2.kelem(i,j,k) = ccomplex_t(eval[0],evec3_t.x);
                                            temp3.kelem(i,j,k) = ccomplex_t(evec3_t.y,evec3_t.z);
                                        }
                                    }//sign
                                } //l3
                            } //l2
                        } //l1
                    } //k
                } //j
            } //i
        }

        D_xx_.kelem(0,0,0) = 1.0;
        D_xy_.kelem(0,0,0) = 0.0;
        D_xz_.kelem(0,0,0) = 0.0;

        D_yy_.kelem(0,0,0) = 1.0;
        D_yz_.kelem(0,0,0) = 0.0;
        D_zz_.kelem(0,0,0) = 0.0;

        //... approximate infinite lattice by inerpolating to sites not convered by current resolution...
        #pragma omp parallel for
        for( size_t i=0; i<D_xx_.size(0); i++ ){
            for( size_t j=0; j<D_xx_.size(1); j++ ){
                for( size_t k=0; k<D_xx_.size(2); k++ ){
                    int ii = (int(i)>nlattice/2)? int(i)-nlattice : int(i);
                    int jj = (int(j)>nlattice/2)? int(j)-nlattice : int(j);
                    int kk = (int(k)>nlattice/2)? int(k)-nlattice : int(k);
                    vec3_t<real_t> kv({real_t(ii),real_t(jj),real_t(kk)});

                    auto align_with_k = [&]( const vec3_t<real_t>& v ) -> vec3_t<real_t>{
                        return v*((v.dot(kv)<0.0)?-1.0:1.0);
                    };

                    vec3_t<real_t> v, l;
                    map_t::iterator it;
                    
                    if( !is_in(i,j,k,mat_invrecip)  ){
                        auto average_lv = [&]( const auto& t1, const auto& t2, const auto& t3, vec3_t<real_t>& v, vec3_t<real_t>& l ) {
                            v = vec3_t<real_t>(0.0); l = vec3_t<real_t>(0.0);
                            int count(0);
                            
                            auto add_lv = [&]( auto it ) -> void {
                                auto q = it->second;++count;
                                l += vec3_t<real_t>({std::real(t1.kelem(q)),std::imag(t1.kelem(q)),std::real(t2.kelem(q))});
                                v += align_with_k(vec3_t<real_t>({std::imag(t2.kelem(q)),std::real(t3.kelem(q)),std::imag(t3.kelem(q))}));
                            };
                            map_t::iterator it;
                            if( (it = iimap.find({ii-1,jj,kk}))!=iimap.end() ){ add_lv(it); }
                            if( (it = iimap.find({ii+1,jj,kk}))!=iimap.end() ){ add_lv(it); }
                            if( (it = iimap.find({ii,jj-1,kk}))!=iimap.end() ){ add_lv(it); }
                            if( (it = iimap.find({ii,jj+1,kk}))!=iimap.end() ){ add_lv(it); }
                            if( (it = iimap.find({ii,jj,kk-1}))!=iimap.end() ){ add_lv(it); }
                            if( (it = iimap.find({ii,jj,kk+1}))!=iimap.end() ){ add_lv(it); }
                            l/=real_t(count); v/=real_t(count);
                        };
                        
                        average_lv(temp1,temp2,temp3,v,l);
                        
                    }else{
                        if( (it = iimap.find({ii,jj,kk}))!=iimap.end() ){
                            auto q = it->second;
                            l = vec3_t<real_t>({std::real(temp1.kelem(q)),std::imag(temp1.kelem(q)),std::real(temp2.kelem(q))});
                            v = align_with_k(vec3_t<real_t>({std::imag(temp2.kelem(q)),std::real(temp3.kelem(q)),std::imag(temp3.kelem(q))}));
                        }
                    }
                    D_xx_.kelem(i,j,k) = l[0];
                    D_xy_.kelem(i,j,k) = l[1];
                    D_xz_.kelem(i,j,k) = l[2];
                    D_yy_.kelem(i,j,k) = v[0];
                    D_yz_.kelem(i,j,k) = v[1];
                    D_zz_.kelem(i,j,k) = v[2];
                }
            }
        }
        
#ifdef PRODUCTION
        #pragma omp parallel for
        for( size_t i=0; i<D_xx_.size(0); i++ ){
            for( size_t j=0; j<D_xx_.size(1); j++ ){
                for( size_t k=0; k<D_xx_.size(2); k++ )
                {
                    vec3_t<real_t> kv = D_xx_.get_k<real_t>(i,j,k);

                    double mu1 = std::real(D_xx_.kelem(i,j,k));
                    // double mu2 = std::real(D_xy_.kelem(i,j,k));
                    // double mu3 = std::real(D_xz_.kelem(i,j,k));

                    vec3_t<real_t> evec1({std::real(D_yy_.kelem(i,j,k)),std::real(D_yz_.kelem(i,j,k)),std::real(D_zz_.kelem(i,j,k))});
                    evec1 /= evec1.norm();

                    // ///////////////////////////////////
                    // // project onto spherical coordinate vectors
                    
                    real_t kr = kv.norm(), kphi = kr>0.0? std::atan2(kv.y,kv.x) : real_t(0.0), ktheta = kr>0.0? std::acos( kv.z / kr ): real_t(0.0);
                    real_t st = std::sin(ktheta), ct = std::cos(ktheta), sp = std::sin(kphi), cp = std::cos(kphi);
                    vec3_t<real_t> e_r( st*cp, st*sp, ct), e_theta( ct*cp, ct*sp, -st), e_phi( -sp, cp, real_t(0.0) );

                    // re-normalise to that longitudinal amplitude is exact
                    double renorm = evec1.dot( e_r ); if( renorm < 0.01 ) renorm = 1.0;

                    // -- store in diagonal components of D_ij
                    D_xx_.kelem(i,j,k) = 1.0;
                    D_yy_.kelem(i,j,k) = evec1.dot( e_theta ) / renorm;
                    D_zz_.kelem(i,j,k) = evec1.dot( e_phi ) / renorm;

                    // spatially dependent correction to vfact = \dot{D_+}/D_+
                    D_xy_.kelem(i,j,k) = 1.0/(0.25*(std::sqrt(1.+24*mu1)-1.));
                }
            }
        }
        D_xy_.kelem(0,0,0) = 1.0;
        D_xx_.kelem(0,0,0) = 1.0;
        D_yy_.kelem(0,0,0) = 0.0;
        D_zz_.kelem(0,0,0) = 0.0;

        // unlink("debug.hdf5");
        // D_xy_.Write_to_HDF5("debug.hdf5","mu1");
        // D_xx_.Write_to_HDF5("debug.hdf5","e1x");
        // D_yy_.Write_to_HDF5("debug.hdf5","e1y");
        // D_zz_.Write_to_HDF5("debug.hdf5","e1z");

#else
        D_xx_.Write_to_HDF5("debug.hdf5","mu1");
        D_xy_.Write_to_HDF5("debug.hdf5","mu2");
        D_xz_.Write_to_HDF5("debug.hdf5","mu3");
        D_yy_.Write_to_HDF5("debug.hdf5","e1x");
        D_yz_.Write_to_HDF5("debug.hdf5","e1y");
        D_zz_.Write_to_HDF5("debug.hdf5","e1z");
#endif   
    }


public:
    // real_t boxlen, size_t ngridother
    explicit lattice_gradient( config_file& the_config, size_t ngridself=64 )
    : boxlen_( the_config.get_value<double>("setup", "BoxLength") ), 
      aini_ ( 1.0/(1.0+the_config.get_value<double>("setup", "zstart")) ),
      ngmapto_( the_config.get_value<size_t>("setup", "GridRes") ), 
      ngrid_( ngridself ), ngrid32_( std::pow(ngrid_, 1.5) ), mapratio_(real_t(ngrid_)/real_t(ngmapto_)),
      //XmL_ ( the_config.get_value<double>("cosmology", "Omega_L") / the_config.get_value<double>("cosmology", "Omega_m") ),
      D_xx_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}), D_xy_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
      D_xz_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}), D_yy_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
      D_yz_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}), D_zz_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
      grad_x_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}), grad_y_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0}),
      grad_z_({ngrid_, ngrid_, ngrid_}, {1.0,1.0,1.0})
    { 
        music::ilog << music::HRULE << std::endl;
        std::string lattice_str = the_config.get_value_safe<std::string>("setup","ParticleLoad","sc");
        const lattice lattice_type = 
            ((lattice_str=="bcc")? lattice_bcc 
            : ((lattice_str=="fcc")? lattice_fcc 
            : ((lattice_str=="rsc")? lattice_rsc 
            : lattice_sc)));

        music::ilog << "PLT corrections for " << lattice_str << " lattice will be computed on " << ngrid_ << "**3 mesh" << std::endl;

        double wtime = get_wtime();
        music::ilog << std::setw(40) << std::setfill('.') << std::left << "Computing PLT eigenmodes "<< std::flush;
        
        init_D( lattice_type );
        // init_D__old();

        music::ilog << std::setw(20) << std::setfill(' ') << std::right << "took " << get_wtime()-wtime << "s" << std::endl;
    }

    inline ccomplex_t gradient( const int idim, std::array<size_t,3> ijk ) const
    {
        real_t ix = ijk[0]*mapratio_, iy = ijk[1]*mapratio_, iz = ijk[2]*mapratio_;

        auto kv = D_xx_.get_k<real_t>( ix, iy, iz );
        auto kmod = kv.norm() / mapratio_ / boxlen_;

        // // project onto spherical coordinate vectors
        auto D_r = std::real(D_xx_.get_cic_kspace({ix,iy,iz}));
        auto D_theta = std::real(D_yy_.get_cic_kspace({ix,iy,iz}));
        auto D_phi = std::real(D_zz_.get_cic_kspace({ix,iy,iz}));
        
        real_t kr = kv.norm(), kphi = kr>0.0? std::atan2(kv.y,kv.x) : 0.0, ktheta = kr>0.0? std::acos( kv.z / kr ) : 0.0;
        real_t st = std::sin(ktheta), ct = std::cos(ktheta), sp = std::sin(kphi), cp = std::cos(kphi);
        
        if( idim == 0 ){
            return ccomplex_t(0.0, kmod*(D_r * st * cp + D_theta * ct * cp - D_phi * sp)); 
        }
        else if( idim == 1 ){
            return ccomplex_t(0.0, kmod*(D_r  * st * sp + D_theta * ct * sp + D_phi * cp)); 
        }
        return ccomplex_t(0.0, kmod*(D_r  * ct - D_theta * st)); 
    }

    inline real_t vfac_corr( std::array<size_t,3> ijk  ) const
    {
        real_t ix = ijk[0]*mapratio_, iy = ijk[1]*mapratio_, iz = ijk[2]*mapratio_;
        const real_t alpha = 1.0/std::real(D_xy_.get_cic_kspace({ix,iy,iz}));
        return 1.0/alpha;
        // // below is for LCDM, but it is a tiny correction for typical starting redshifts:
        //! X = \Omega_\Lambda / \Omega_m
        // return 1.0 / (alpha - (2*std::pow(aini_,3)*alpha*(2 + alpha)*XmL_*Hypergeometric2F1((3 + alpha)/3.,(5 + alpha)/3.,
        //     (13 + 4*alpha)/6.,-(std::pow(aini_,3)*XmL_)))/
        //     ((7 + 4*alpha)*Hypergeometric2F1(alpha/3.,(2 + alpha)/3.,(7 + 4*alpha)/6.,-(std::pow(aini_,3)*XmL_))));
    }

};

}
