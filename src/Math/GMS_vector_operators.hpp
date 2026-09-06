
/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef GMS_VECTOR_OPERATORS_H_
#define GMS_VECTOR_OPERATORS_H_


#include <valarray>
#include <vector>
#include <functional>
#include <algorithm>

namespace   gms {

	// TODO: Add SIMD vectorization.

	


		/*

		           ***   Declarations of namespace operators std::vector and std::valarray wrapper operators.  ***
		*/

		    
		/*
		               ***  std::vector tag declaration/definition.
		*/

		                      struct   tag_std_vector;

							  struct   tag_std_valarray;

							  /*
							         Primary template
							  */
							  template<class Vectorable> struct tag;
							  
							  /*
							          std::vector<T,Alloc> specialization
							  */
							  template<> struct tag<std::vector<float>>{

								  typedef std::vector<float> v_type;
							  };

							  /*
							         Fills Vector Out with the result of call to Fx
							  */
							  template<typename Vector, class Fx> struct
								  data_filler{

								  __forceinline  std::enable_if<std::is_same<Vector,Vector>::value, void> operator()(Vector const &In, Vector &&Out, Fx fx) {
									  static_assert(In.size() == Out.size(), "Invalid size: In.size() != Out.size()");
									  for (std::size_t i{ 0 }; i != In.size(); ++i) {

										  Out.operator[](i) = fx(In.operator[](i));
									  }
								  }
								  

							  };
								  
							  

							  /*
							      
							  */
							  template<> struct tag<std::valarray<float>>{

								  typedef tag_std_valarray v_type;
							  };

							          /*

									   Creates operator/ for dividing  vector by the  scalar of the same type.
									   No attempt for division by the zero check is made.
									   */
				template<typename T> auto
				operator /(std::vector<T> const &v, T const s)->std::vector<decltype(T{} / T{})>{

				std::vector<T> v_res(v.size());
				std::divides<T> divides;
				for (std::size_t i{ 0 }; i != v_res.size(); ++i){

					v_res.operator[](i) = divides.operator()(v.operator[](i), s);
				}

				return v_res;
			}


				

			/*

			          Creates an operator/ for dividing std::vector<T1> by the scalar T2. 
			          No attempt for division by the zero check is made.                         
			  
			*/
			template<typename T1, typename T2> typename std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type 
				operator / (std::vector<T1> const &v, T2 const s){

					static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
					std::vector<T1> v_res(v.size());
					std::divides<T1> divides;
					for (std::size_t i{ 0 }; i != v_res.size(); ++i){

						v_res.operator[](i) = divides.operator()(v.operator[](i), static_cast<T1>(s));
					}
					return v_res;
				}

			/*

			Creates operator/ for dividing  vector by the  scalar of the same type(swapped arguments).
			No attempt for division by the zero check is made.

			*/
			template<typename T>  auto
				operator/(T const s, std::vector<T> const &v)->std::vector<decltype(T{} / T{})>{

					std::vector<T> v_res(v.size());
					std::divides<T> divides;
					for (std::size_t i{ 0 }; i != v_res.size(); ++i){

						v_res.operator[](i) = divides.operator()(v.operator[](i), s);
					}
					return v_res;
				}

				/*

				   Creates an operator/ for dividing std::vector<T1> by the scalar T2(swapped arguments). 
			          No attempt for division by the zero check is made. 

				*/
				template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
					operator/(T2 const s, std::vector<T1> const &v){

						static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
						std::vector<T1> v_res(v.size());
						std::divides<T1> divides;
						for (std::size_t i{ 0 }; i != v_res.size(); ++i){

							v_res.operator[](i) = divides.operator()(v.operator[](i), static_cast<T1>(s));
						}
						return v_res;
					}

				/*

						Creates an operator/ for dividing std::vector<T> by the vector<T>.
						No attempt for division by the zero check is made.

						*/
				template<typename T> auto
					operator/(std::vector<T> const &v1, std::vector<T> const &v2)->std::vector<decltype(T{} / T{})> {

						BOOST_ASSERT_MSG(v1.size() == v2.size(), "v1.size() != v2.size()");
						std::vector<T> v_res(v1.size());
						std::divides<T> divides;
						for (std::size_t i{ 0 }; i != v_res.size(); ++i){

							v_res.operator[](i) = divides.operator()(v1.operator[](i), v2.operator[](i));
						}
						return v_res;
					}

					/*

					      Creates an operator/ for dividing std::vector<T1> by the vector<T2>. 
			            No attempt for division by the zero check is made.   

					*/
				template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
					operator/(std::vector<T1> const &v1, std::vector<T2> const &v2){

						static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
						BOOST_ASSERT_MSG(v1.size() == v2.size(), "v1.size() != v2.size()");
						std::vector<T1> v_res(v1.size());
						std::divides<T1> divides;
						for (std::size_t i{ 0 }; i != v_res.size(); ++i){

							v_res.operator[](i) = divides.operator()(v1.operator[](i), static_cast<T1>(v2.operator[](i)));
						}
						return v_res;
					}


				/*

						Creates operator* for multiplying vector<T> by the  scalar T.


						*/
				template<typename T>   auto
					operator*(std::vector<T> const &v, T const s)->std::vector<decltype(T{}*T{})>{

						std::vector<T> v_res(v.size());
						std::multiplies<T> mul;
						for (std::size_t i{ 0 }; i != v_res.size(); ++i){

							v_res.operator[](i) = mul.operator()(v.operator[](i), s);
						}
						return v_res;
					}

					/*

					Creates an operator* for multiplying std::vector<T1> by the scalar T2.
					

					*/
					template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
						operator*(std::vector<T1> const &v, T2 const s){

							static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
							std::vector<T1> v_res(v.size());
							std::multiplies<T1> mul;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i){

								v_res.operator[](i) = mul.operator()(v.operator[])(i), static_cast<T1>(s));
							}
							return v_res;
						}

					/*

					Creates operator* for multiplying  sclar T by the  vector of the same type(swapped arguments).
					

					*/

					template<typename T>  auto
						operator*(T const s, std::vector<T> const &v)->std::vector<decltype(T{} *T{})> {

							std::vector<T> v_res(v.size());
							std::multiplies<T> mul;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = mul.operator()(v.operator[](i), s);
							}
							return v_res;
						}

					/*

					Creates an operator* for multiplying std::vector<T1> by the scalar T2(swapped arguments).
					

					*/

					template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
						operator*(T2 const s, std::vector<T1> const &v){

							static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
							std::vector<T1> v_res(v.size());
							std::multiplies<T1> mul;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = mul.operator()(v.operator[](i), static_cast<T1>(s));
							}
							return v_res;
						}

					/*

					Creates an operator* for multiplying std::vector<T> by the vector<T>.
					

					*/
					template<typename T>  auto
						operator*(_In_ std::vector<T> const &v1, _In_ std::vector<T> const &v2)->std::vector<decltype(T{} *T{})>{

							BOOST_ASSERT_MSG(v1.size() == v2.size(), "Invalid size: v1.size() != v2.size()");
							std::vector<T> v_res(v1.size());
							std::multiplies<T> mul;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = mul.operator()(v1.operator[](i), v2.operator[](i));
							}
							return v_res;
						}

					/*

					Creates an operator* for multiplying std::vector<T1> by the vector<T2>.
					

					*/
					template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
						operator*(_In_ std::vector<T1> const &v1, _In_ std::vector<T2> const &v2) {

							static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
							BOOST_ASSERT_MSG(v1.size() == v2.size(), "Invalid size: v1.size() != v2.size()");
							std::vector<T1> v_res(v1.size());
							std::multiplies<T1> mul;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = mul.operator()(v1.operator[](i), static_cast<T1>(v2.operator[](i)));
							}
							return v_res;

						}

					/*

					Creates operator+ for adding vector<T> by the  scalar T.


					*/
					template<typename T> auto
						operator+(_In_ std::vector<T> const &v, T const s)->std::vector<decltype(T{}+T{})> {

							std::vector<T> v_res(v.size());
							std::plus<T> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v.operator[](i), s);
							}
							return v_res;
						}

					/*

					Creates an operator+ for addidtion of std::vector<T1> by the scalar T2.


					*/

					template<typename T1, typename T2> std::enable_if<!std::is_same<T1, T2>::value, std::vector<T1>>::type
						operator+(_In_ std::vector<T1> const &v, T2 const s) {
							static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
							std::vector<T1> v_res(v.size());
							std::plus<T1> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v.operator[](i), static_cast<T1>(s));
							}
						
							return v_res;
						}

					/*

					Creates operator+ for addition  scalar T by the  vector of the same type(swapped arguments).


					*/
					template<typename T> auto
						operator+(_In_ T const s, _In_ std::vector<T> const &v)->std::vector < decltype(T{} +T{})> {

							std::vector<T> v_res(v.size());
							std::plus<T> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v.operator[](i), s);
							}
							return v_res;
						}
					
					/*

					Creates an operator* for addition of std::vector<T1> by the scalar T2(swapped arguments).


					*/
					template<typename T1, typename T2> std::enable_if<!std::is_same<T1,T2>::value,std::vector<T1>>::type
						operator+(_In_ T2 const s, _In_ std::vector<T1> const &v) {

							static_assert(std::is_convertible<T2, T1>::operator bool, "Cannot convert from T2 to T1");
							std::vector<T1> v_res(v.size());
							std::plus<T1> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v.operator[](i), static_cast<T2>(s));
							}
							return v_res;
						}

					/*

					Creates an operator+ for adding std::vector<T> by the vector<T>.


					*/
					template<typename T> auto
						operator+(_In_ std::vector<T> const &v1, _In_ std::vector<T> const &v2)->std::vector<decltype(T{}+T{})> {

							BOOST_ASSERT_MSG(v1.size() == v2.size(), "Invalid size: v1.size() != v2.size()");
							std::vector<T> v_res(v1.size());
							std::plus<T> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v1.operator[](i), v2.operator[](i));
							}
							return v_res;
						}

					/*

					Creates an operator+ for addition of  std::vector<T1> by the vector<T2>.


					*/
					template<typename T1, typename T2> std::enable_if<!(std::is_same<T1,T2>::value) && (std::is_convertible<T2,T1>::value),std::vector<T1>>::type
						operator+(_In_ std::vector<T1> const &v1, _In_ std::vector<T2> const v2) {
							
							BOOST_ASSERT_MSG(v1.size() == v2.size(), "Invalid size: v1.size() != v2.size()");
							std::vector<T1> v_res(v1.size());
							std::plus<T1> add;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = add.operator()(v1.operator[](i), static_cast<T1>(v2.operator[](i)));
							}
							return v_res;
						}

					/*

					Creates operator- for subtracting vector<T> by the  scalar T.


					*/
					template<typename T> auto
						operator-(_In_ std::vector<T> const &v, _In_ T const s)->std::vector < decltype(T{}-T{})> {

							std::vector<T> v_res(v.size());
							std::minus<T> sub;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = sub.operator()(v_res.operator[](i), s);
							}
							return v_res;
						}

					/*

					Creates operator- for subtracting vector<T1> by the  scalar T2.


					*/
					template<typename T1, typename T2> std::enable_if<!(std::is_same<T1,T2>::value) && (std::is_convertible<T2,T1>::value),std::vector<T1>>::type
						operator-(_In_ std::vector<T1> const &v, _In_ T const s) {

							std::vector<T1> v_res(v.size());
							std::minus<T1> sub;
							for (std::size_t i{ 0 }; i != v_res.size(); ++i) {

								v_res.operator[](i) = sub.operator()(v.operator[](i), static_cast<T1>(s));
							}
							return v_res;
						}
						
			
	

	

		                   
		/*
		               SIMD AVX  version of addVector function
					   Primary template.
		*/
		template<typename T> void 
			addAVXVector(_Out_ std::vector<T> &, _In_ std::vector<T> const &, std::vector<T> const &);

		

		/*
		                Float specialization. Do not check for c.size() % 8 == 0 || b.size() % 8 == 0 || a.size() % 8 == 0
						Non-temporal streams in use.
		*/
		template<>  void  
			addAVXVector<float>(_Out_ std::vector<float> &c, _In_ std::vector<float> const &b, _In_ std::vector<float> const &a){

				for (std::size_t i_idx{ 0 }; i_idx != c.size(); i_idx += 8){

					_mm256_stream_ps(&c.operator[](i_idx), _mm256_add_ps(_mm256_loadu_ps(&b.operator[](i_idx)), _mm256_loadu_ps(&a.operator[](i_idx))));
				}
			}

		


		

		/*

		               Double specialization. Do not check for c.size() % 4 == 0 || b.size() % 4 == 0 ||  a.size() % 4 == 0
		*/
		template<>  void
			addAVXVector<double>(_Out_ std::vector<double> &c, _In_ std::vector<double> const &b, _In_ std::vector<double> const &a) {

				for (std::size_t i_idx{ 0 }; i_idx != c.size(); i_idx += 4) {

					_mm256_stream_pd(&c.operator[](i_idx), _mm256_add_pd(_mm256_loadu_pd(&b.operator[](i_idx)), _mm256_loadu_pd(&a.operator[](i_idx))));
				}
			}
		  
		/*

		               SIMD SSE  version of addVector function
					   Primary template.

		*/
		template<typename T>  void
			addSSEVector(_Out_ std::vector<T> &, _In_ std::vector<T> const &, _In_ std::vector<T> const &);

		/*

		               Float specialization. Do not check for c.size() % 4 == 0 || b.size() % 4 == 0 ||  a.size() % 4 == 0

		*/
		template<>    void
			addSSEVector<float>(_Out_ std::vector<float> &c, _In_ std::vector<float> const &b, _In_ std::vector<float> const &a) {

				for (std::size_t i_idx{ 0 }; i_idx != c.size(); i_idx += 4){

					_mm_stream_ps(&c.operator[](i_idx), _mm_add_ps(_mm_loadu_ps(&b.operator[](i_idx)), _mm_loadu_ps(&a.operator[](i_idx))));
				}
			}

		/*

		               Double specialization. Do not check for c.size() % 2 == 0 || b.size() % 2 == 0 ||  a.size() % 2 == 0.
					   Non-temporal stream in use.

		*/
		template<>    void
			addSSEVector<double>(_Out_ std::vector<double> &c, _In_ std::vector<double> const &b, _In_ std::vector<double> const &a) {

				for (std::size_t i_idx{ 0 }; i_idx != c.size(); i_idx += 2) {

					_mm_stream_pd(&c.operator[](i_idx), _mm_add_pd(_mm_loadu_pd(&b.operator[](i_idx)), _mm_loadu_pd(&a.operator[](i_idx))));
				}
			}

		/*
		               std::vector<T> a,b,c , where   c = a + b. No check is made on the vectors length.
					   Should be faster than equivalent operator+ based version.
		*/
		    template<typename T>   void    
				addVector(_Out_ std::vector<T>  &c, _In_ std::vector<T> const &b, _In_ std::vector<T> const &a){

					for (std::size_t i_idx{ 0 }; i_idx != a.size(); ++i_idx) {

						c.operator[](i_idx) = b.operator[](i_idx) + a.operator[](i_idx);
					}
		}

			/*
			           std::vector<T> a,b, scalar: s where a = b + s. No check is made on the vectors length.
					   Should be faster than equivalent operator+ based version.
			
			*/
			template<typename T>   void
				addVector(_Out_ std::vector<T> &c, _In_ std::vector<T> const &b, _In_ T const s) {

					for (std::size_t i_idx{ 0 }; i_idx != b.size(); ++i_idx) {

						c.operator[](i_idx) = b.operator[](i) + s;
					}
				}

			/*
			std::vector<T> a,b, scalar: s where a = b + s(swapped). No check is made on the vectors length.
			Should be faster than equivalent operator+ based version.

			*/
			template<typename T>   void
				addVector(_Out_ std::vector<T> &c, _In_ T const s, _In_ std::vector<T> const &b) {

					for (std::size_t i_idx{ 0 }; i_idx != b.size(); ++i_idx) {

						c.operator[](i_idx) = b.operator[](i_idx) + s;
					}
				}
	
}
#endif /*GMS_VECTOR_OPERATORS_H_*/
