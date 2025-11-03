/* Copyright (c) 2015, Bernard Gingold. License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Generic templated Matrix 4x4  classes based upon SIMD types __m256d and __m256 -  classes implementation.


*/

const __m128 gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::UTMask[4] = 
{ _mm_set_ps(1.f,1.f,1.f,1.f), _mm_set_ps(0.f,1.f,1.f,1.f), _mm_set_ps(0.f,0.f,1.f,1.f), _mm_set_ps(0.f,0.f,0.f,1.f) };

const __m128 gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::LTMask[4] =
{ _mm_set_ps(1.f, 0.f, 0.f, 0.f), _mm_set_ps(1.f, 1.f, 0.f, 0.f), _mm_set_ps(1.f, 1.f, 1.f, 0.f), _mm_set_ps(1.f, 1.f, 1.f, 1.f) };

 gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::SIMDMat4x4(){

	this->m_fMat4x4[0] = _mm_setzero_ps();
	this->m_fMat4x4[1] = _mm_setzero_ps();
	this->m_fMat4x4[2] = _mm_setzero_ps();
	this->m_fMat4x4[3] = _mm_setzero_ps();
}



 gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::SIMDMat4x4(_In_ const ArgsHolder<float> &args){

	this->m_fMat4x4[0] = _mm_set_ps(args.a00, args.a01, args.a02, args.a03);
	this->m_fMat4x4[1] = _mm_set_ps(args.a10, args.a11, args.a12, args.a13);
	this->m_fMat4x4[2] = _mm_set_ps(args.a20, args.a21, args.a22, args.a23);
	this->m_fMat4x4[3] = _mm_set_ps(args.a30, args.a31, args.a32, args.a33);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::SIMDMat4x4(_In_ const __m128 &v){

	this->m_fMat4x4[0] = v;
	this->m_fMat4x4[1] = v;
	this->m_fMat4x4[2] = v;
	this->m_fMat4x4[3] = v;
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::SIMDMat4x4(_In_ std::initializer_list<__m128> const &arg)
{
	this->m_fMat4x4[0] = arg.begin()[0];
	this->m_fMat4x4[1] = arg.begin()[1];
	this->m_fMat4x4[2] = arg.begin()[2];
	this->m_fMat4x4[3] = arg.begin()[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::SIMDMat4x4(_In_ const float s){

	this->m_fMat4x4[0] = _mm_set1_ps(s);
	this->m_fMat4x4[1] = _mm_set1_ps(s);
	this->m_fMat4x4[2] = _mm_set1_ps(s);
	this->m_fMat4x4[3] = _mm_set1_ps(s);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::SIMDMat4x4(_In_ float(&ar)[16]){

	constexpr int N{ 4 };

	for (int i{ 0 }; i != N; ++i){
		for (int j{ 0 }; j != N; ++j){
			this->m_fMat4x4[i].m128_f32[j] = ar[i + N * j];
		}
	}
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::SIMDMat4x4(_In_ __m128(&var)[4]){

	this->m_fMat4x4[0] = var[0];
	this->m_fMat4x4[1] = var[1];
	this->m_fMat4x4[2] = var[2];
	this->m_fMat4x4[3] = var[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::SIMDMat4x4(_In_ SIMDMat4x4 const &rhs) {

	this->m_fMat4x4[0] = rhs.m_fMat4x4[0];
	this->m_fMat4x4[1] = rhs.m_fMat4x4[1];
	this->m_fMat4x4[2] = rhs.m_fMat4x4[2];
	this->m_fMat4x4[3] = rhs.m_fMat4x4[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::SIMDMat4x4(_In_ SIMDMat4x4 &&rhs){

	this->m_fMat4x4[0] = std::move(rhs.m_fMat4x4[0]);
	this->m_fMat4x4[1] = std::move(rhs.m_fMat4x4[1]);
	this->m_fMat4x4[2] = std::move(rhs.m_fMat4x4[2]);
	this->m_fMat4x4[3] = std::move(rhs.m_fMat4x4[3]);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &   
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator=(_In_ SIMDMat4x4 const &rhs){

	if (this == &rhs) return (*this);

	this->m_fMat4x4[0] = rhs.m_fMat4x4[0];
	this->m_fMat4x4[1] = rhs.m_fMat4x4[1];
	this->m_fMat4x4[2] = rhs.m_fMat4x4[2];
	this->m_fMat4x4[3] = rhs.m_fMat4x4[3];
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &    
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator=(_In_ SIMDMat4x4 &&rhs){

	if (this == &rhs) return (*this);

	this->m_fMat4x4[0] = std::move(rhs.m_fMat4x4[0]);
	this->m_fMat4x4[1] = std::move(rhs.m_fMat4x4[1]);
	this->m_fMat4x4[2] = std::move(rhs.m_fMat4x4[2]);
	this->m_fMat4x4[3] = std::move(rhs.m_fMat4x4[3]);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &    
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator+=(_In_ SIMDMat4x4 const &rhs){

	this->m_fMat4x4[0] = _mm_add_ps(this->m_fMat4x4[0], rhs.m_fMat4x4[0]);
	this->m_fMat4x4[1] = _mm_add_ps(this->m_fMat4x4[1], rhs.m_fMat4x4[1]);
	this->m_fMat4x4[2] = _mm_add_ps(this->m_fMat4x4[2], rhs.m_fMat4x4[2]);
	this->m_fMat4x4[3] = _mm_add_ps(this->m_fMat4x4[3], rhs.m_fMat4x4[3]);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &    
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator+=(_In_ __m128 const &v){

	this->m_fMat4x4[0] = _mm_add_ps(this->m_fMat4x4[0], v);
	this->m_fMat4x4[1] = _mm_add_ps(this->m_fMat4x4[1], v);
	this->m_fMat4x4[2] = _mm_add_ps(this->m_fMat4x4[2], v);
	this->m_fMat4x4[3] = _mm_add_ps(this->m_fMat4x4[3], v);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator+=(_In_ const float s){

	this->m_fMat4x4[0] = _mm_add_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	this->m_fMat4x4[1] = _mm_add_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	this->m_fMat4x4[2] = _mm_add_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	this->m_fMat4x4[3] = _mm_add_ps(this->m_fMat4x4[3], _mm_set1_ps(s));
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator-=(_In_ SIMDMat4x4 const &rhs){

	this->m_fMat4x4[0] = _mm_sub_ps(this->m_fMat4x4[0], rhs.m_fMat4x4[0]);
	this->m_fMat4x4[1] = _mm_sub_ps(this->m_fMat4x4[1], rhs.m_fMat4x4[1]);
	this->m_fMat4x4[2] = _mm_sub_ps(this->m_fMat4x4[2], rhs.m_fMat4x4[2]);
	this->m_fMat4x4[3] = _mm_sub_ps(this->m_fMat4x4[3], rhs.m_fMat4x4[3]);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator-=(_In_ __m128 const &v){

	this->m_fMat4x4[0] = _mm_sub_ps(this->m_fMat4x4[0], v);
	this->m_fMat4x4[1] = _mm_sub_ps(this->m_fMat4x4[1], v);
	this->m_fMat4x4[2] = _mm_sub_ps(this->m_fMat4x4[2], v);
	this->m_fMat4x4[3] = _mm_sub_ps(this->m_fMat4x4[3], v);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator-=(_In_ const float s){

	this->m_fMat4x4[0] = _mm_sub_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	this->m_fMat4x4[1] = _mm_sub_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	this->m_fMat4x4[2] = _mm_sub_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	this->m_fMat4x4[3] = _mm_sub_ps(this->m_fMat4x4[3], _mm_set1_ps(s));
	return (*this);
	
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator*=(_In_ SIMDMat4x4  &rhs){

	rhs.transpose();
	__m128 t0(_mm_mul_ps(this->m_fMat4x4[0], rhs.m_fMat4x4[0]));
	this->m_fMat4x4[0] = _mm_add_ps(t0, _mm_permute_ps(t0, 21));
	__m128 t1(_mm_mul_ps(this->m_fMat4x4[1], rhs.m_fMat4x4[1]));
	this->m_fMat4x4[1] = _mm_add_ps(t1, _mm_permute_ps(t1, 21));
	__m128 t2(_mm_mul_ps(this->m_fMat4x4[2], rhs.m_fMat4x4[2]));
	this->m_fMat4x4[2] = _mm_add_ps(t2, _mm_permute_ps(t2, 21));
	__m128 t3(_mm_mul_ps(this->m_fMat4x4[3], rhs.m_fMat4x4[3]));
	this->m_fMat4x4[3] = _mm_add_ps(t3, _mm_permute_ps(t3, 21));
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &      
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator*=(_In_ __m128 const &v){

	this->m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], v);
	this->m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], v);
	this->m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], v);
	this->m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], v);
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &       
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator*=(_In_ const float s){

	this->m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	this->m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	this->m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	this->m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], _mm_set1_ps(s));
	return (*this);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float> &        
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator/=(_In_ const float s){

	this->m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	this->m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	this->m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	this->m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], _mm_set1_ps(s));
	
	return (*this);
}

__m128                    gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator[](_In_ const int index){
	constexpr int SIZE{ 3 };
	_ASSERT((index >= 0) && (index <= SIZE));
	return this->m_fMat4x4[index];
}

const  __m128             gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator[](_In_ const int index) const{

	constexpr int SIZE{ 3 };
	_ASSERT((index >= 0) && (index <= SIZE));
	return this->m_fMat4x4[index];
}

std::ostream &        gms::math::operator<<(_In_ std::ostream &os, const SIMDMat4x4<gms::math::ArgsHolder<float>, float> &rhs) {

	constexpr int const N{ 4 };
	for (int i{ 0 }; i != N; ++i){
		for (int j{ 0 }; j != N; ++j){
			os << rhs.m_fMat4x4[i].m128_f32[j] << std::endl;
		}
		os << std::endl;
	}
	return os;
}


gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator==(_In_ SIMDMat4x4 const &rhs){

	
	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpeq_ps(this->m_fMat4x4[0], rhs.m_fMat4x4[0]);
	ret.m_fMat4x4[1] = _mm_cmpeq_ps(this->m_fMat4x4[1], rhs.m_fMat4x4[1]);
	ret.m_fMat4x4[2] = _mm_cmpeq_ps(this->m_fMat4x4[2], rhs.m_fMat4x4[2]);
	ret.m_fMat4x4[3] = _mm_cmpeq_ps(this->m_fMat4x4[3], rhs.m_fMat4x4[3]);

	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator==(_In_ __m128 const &rhs){

	
	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpeq_ps(this->m_fMat4x4[0], rhs);
	ret.m_fMat4x4[1] = _mm_cmpeq_ps(this->m_fMat4x4[1], rhs);
	ret.m_fMat4x4[0] = _mm_cmpeq_ps(this->m_fMat4x4[2], rhs);
	ret.m_fMat4x4[0] = _mm_cmpeq_ps(this->m_fMat4x4[3], rhs);

	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator==(_In_ const float s){

	
	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpeq_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	ret.m_fMat4x4[1] = _mm_cmpeq_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	ret.m_fMat4x4[2] = _mm_cmpeq_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	ret.m_fMat4x4[3] = _mm_cmpeq_ps(this->m_fMat4x4[3], _mm_set1_ps(s));

	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator!=(_In_ SIMDMat4x4 const &rhs){


	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpneq_ps(this->m_fMat4x4[0], rhs.m_fMat4x4[0]);
	ret.m_fMat4x4[1] = _mm_cmpneq_ps(this->m_fMat4x4[1], rhs.m_fMat4x4[1]);
	ret.m_fMat4x4[2] = _mm_cmpneq_ps(this->m_fMat4x4[2], rhs.m_fMat4x4[2]);
	ret.m_fMat4x4[3] = _mm_cmpneq_ps(this->m_fMat4x4[3], rhs.m_fMat4x4[3]);

	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator!=(_In_ __m128 const &rhs){

	
	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpneq_ps(this->m_fMat4x4[0], rhs);
	ret.m_fMat4x4[1] = _mm_cmpneq_ps(this->m_fMat4x4[1], rhs);
	ret.m_fMat4x4[2] = _mm_cmpneq_ps(this->m_fMat4x4[2], rhs);
	ret.m_fMat4x4[3] = _mm_cmpneq_ps(this->m_fMat4x4[3], rhs);

	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>         
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::operator!=(_In_ const float s){

	SIMDMat4x4 ret;
	ret.m_fMat4x4[0] = _mm_cmpneq_ps(this->m_fMat4x4[0], _mm_set1_ps(s));
	ret.m_fMat4x4[1] = _mm_cmpneq_ps(this->m_fMat4x4[1], _mm_set1_ps(s));
	ret.m_fMat4x4[2] = _mm_cmpneq_ps(this->m_fMat4x4[2], _mm_set1_ps(s));
	ret.m_fMat4x4[3] = _mm_cmpneq_ps(this->m_fMat4x4[3], _mm_set1_ps(s));

	return (ret);
}

void                    gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::insert_at(_In_  int index, _In_ __m128 &elem) {

	_ASSERT((index >= 0) && (index <= 3));
	this->m_fMat4x4[index] = elem;
	
}

void                    gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::transpose(){

	// Perform Matrix transposition in-place.
	_MM_TRANSPOSE4_PS(this->m_fMat4x4[0], this->m_fMat4x4[1], this->m_fMat4x4[2], this->m_fMat4x4[3]);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>           gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::transpose() const{

	this->transpose();
	SIMDMat4x4 mat4x4;
	mat4x4.m_fMat4x4[0] = this->m_fMat4x4[0];
	mat4x4.m_fMat4x4[1] = this->m_fMat4x4[1];
	mat4x4.m_fMat4x4[2] = this->m_fMat4x4[2];
	mat4x4.m_fMat4x4[3] = this->m_fMat4x4[3];
	return (mat4x4);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>           gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>::identity(){

	SIMDMat4x4 mat4x4;
	mat4x4.m_fMat4x4[0] = _mm_set_ps(0.f, 0.f, 0.f, 1.f);
	mat4x4.m_fMat4x4[1] = _mm_set_ps(0.f, 0.f, 1.f, 0.f);
	mat4x4.m_fMat4x4[2] = _mm_set_ps(0.f, 1.f, 0.f, 0.f);
	mat4x4.m_fMat4x4[3] = _mm_set_ps(1.f, 0.f, 0.f, 0.f);
	
	return (mat4x4);
}

void                       gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::upper_triangular(){
	

	// In-place modification.
	MaskUT<float> UT;
	this->m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], UT.v0);
	this->m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], UT.v1);
	this->m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], UT.v2);
	this->m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], UT.v3);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>        gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::upper_triangular() const{

	MaskUT<float> UT;
	SIMDMat4x4 mat4x4;
	mat4x4.m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], UT.v0);
	mat4x4.m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], UT.v1);
	mat4x4.m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], UT.v2);
	mat4x4.m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], UT.v3);
	return (mat4x4);
}

void                       gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::lower_triangular(){

	// In place modification.
	MaskLT<float> LT;
	this->m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], LT.v0);
	this->m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], LT.v0);
	this->m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], LT.v2);
	this->m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], LT.v3);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>         
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::lower_triangular() const{

	MaskLT<float> LT;
	SIMDMat4x4 mat4x4;
	mat4x4.m_fMat4x4[0] = _mm_mul_ps(this->m_fMat4x4[0], LT.v0);
	mat4x4.m_fMat4x4[1] = _mm_mul_ps(this->m_fMat4x4[1], LT.v1);
	mat4x4.m_fMat4x4[2] = _mm_mul_ps(this->m_fMat4x4[2], LT.v2);
	mat4x4.m_fMat4x4[3] = _mm_mul_ps(this->m_fMat4x4[3], LT.v3);
	return (mat4x4);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>   
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::symmetric_matrix(){

	// Copy *this to temporary.
	SIMDMat4x4 temp = SIMDMat4x4<gms::math::ArgsHolder<float>,float>(*this);
	this->transpose();
	temp.m_fMat4x4[0] = _mm_cmpeq_ps(temp.m_fMat4x4[0], this->m_fMat4x4[0]);
	temp.m_fMat4x4[1] = _mm_cmpeq_ps(temp.m_fMat4x4[1], this->m_fMat4x4[1]);
	temp.m_fMat4x4[2] = _mm_cmpeq_ps(temp.m_fMat4x4[2], this->m_fMat4x4[2]);
	temp.m_fMat4x4[3] = _mm_cmpeq_ps(temp.m_fMat4x4[3], this->m_fMat4x4[3]);
	// return temp modified.
	return temp;
}

float                          gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::trace() const{

	float f_trace = this->m_fMat4x4[0].m128_f32[0] + this->m_fMat4x4[1].m128_f32[1] +
		this->m_fMat4x4[2].m128_f32[2] + this->m_fMat4x4[3].m128_f32[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::x_rotation(_In_ const float f_rad){

	// Create X rotation Matrix.
	float sin_theta{ ::sin(f_rad) };
	float cos_theta{ ::cos(f_rad) };
	ArgsHolder<float> ah;
	ah.a00 = 1.f, ah.a01 = 0x00000000, ah.a02 = 0x00000000, ah.a03 = 0x00000000;
	ah.a10 = 0x00000000, ah.a11 = cos_theta, ah.a12 = sin_theta, ah.a13 = 0x00000000;
	ah.a20 = 0x00000000, ah.a21 = -sin_theta, ah.a22 = cos_theta, ah.a23 = 0x00000000;
	ah.a30 = 0x00000000, ah.a31 = 0x00000000, ah.a32 = 0x00000000, ah.a33 = 1.f;
	SIMDMat4x4 ret = SIMDMat4x4{ ah };
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::y_rotation(_In_ const float f_rad){

	// Create Y rotation Matrix.
	float sin_theta{ ::sin(f_rad) };
	float cos_theta{ ::cos(f_rad) };
	ArgsHolder<float> ah;
	ah.a00 = cos_theta, ah.a01 = 0x00000000, ah.a02 - sin_theta, ah.a03 = 0x00000000;
	ah.a10 = 0x00000000, ah.a11 = 1.f, ah.a12 = 0x00000000, ah.a13 = 0x00000000;
	ah.a20 = sin_theta, ah.a21 = 0x00000000, ah.a22 = cos_theta, ah.a23 = 0x00000000;
	ah.a30 = 0x00000000, ah.a31 = 0x00000000, ah.a32 = 0x00000000, ah.a33 = 1.f;
	SIMDMat4x4 ret = SIMDMat4x4{ ah };
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>     
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::z_rotation(_In_ const float f_rad){

	// Create Z rotation Matrix.
	float sin_theta{ ::sin(f_rad) };
	float cos_theta{ ::cos(f_rad) };
	ArgsHolder<float> ah;
	ah.a00 = cos_theta, ah.a01 = sin_theta, ah.a02 = 0x00000000, ah.a03 = 0x00000000;
	ah.a10 = -sin_theta, ah.a11 = cos_theta, ah.a12 = 0x00000000, ah.a13 = 0x00000000;
	ah.a20 = 0x00000000, ah.a21 = 0x00000000, ah.a22 = 1.f, ah.a23 = 0x00000000;
	ah.a30 = 0x00000000, ah.a31 = 0x00000000, ah.a32 = 0x00000000, ah.a33 = 1.f;
	SIMDMat4x4 ret = SIMDMat4x4{ ah };
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>      
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::scale(_In_ const float f_s){

	ArgsHolder<float> ah;
	ah.a00 = f_s, ah.a01 = 0x00000000, ah.a02 = 0x00000000, ah.a03 = 0x00000000;
	ah.a10 = 0x00000000, ah.a11 = f_s, ah.a12 = 0x00000000, ah.a13 = 0x00000000;
	ah.a20 = 0x00000000, ah.a21 = 0x00000000, ah.a22 = f_s, ah.a23 = 0x00000000;
	ah.a30 = 0x00000000, ah.a31 = 0x00000000, ah.a32 = 0x00000000, ah.a33 = 1.f;//1.f
	SIMDMat4x4 ret = SIMDMat4x4{ ah };
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>      
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::translate(_In_ const float f_x, _In_ const float f_y, _In_ const float f_z){

	ArgsHolder<float> ah;
	ah.a00 = 1.f, ah.a01 = 0x00000000, ah.a02 = 0x00000000, ah.a03 = 0x00000000;
	ah.a10 = 0x00000000, ah.a11 = 1.f, ah.a12 = 0x00000000, ah.a13 = 0x00000000;
	ah.a20 = 0x00000000, ah.a21 = 0x00000000, ah.a22 = 1.f, ah.a23 = 0x00000000;
	ah.a30 = f_x, ah.a31 = f_y, ah.a32 = f_z, ah.a33 = 1.f;
	
	SIMDMat4x4 ret = SIMDMat4x4{ ah };
	return (ret);
}

__m128     gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::row_0() const{

	return (this->m_fMat4x4[0]);
}

__m128     gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::row_1() const{

	return (this->m_fMat4x4[1]);
}

__m128     gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::row_2() const{

	return (this->m_fMat4x4[2]);
}

__m128     gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>,float>::row_3() const{

	return (this->m_fMat4x4[3]);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>         
gms::math::operator+(_In_ const SIMDMat4x4<gms::math::ArgsHolder<float>, float> &lhs, _In_ const SIMDMat4x4<gms::math::ArgsHolder<float>,float> &rhs){

	SIMDMat4x4<gms::math::ArgsHolder<float>,float> ret;
	ret.m_fMat4x4[0] = _mm_add_ps(lhs.m_fMat4x4[0], rhs.m_fMat4x4[0]);
	ret.m_fMat4x4[1] = _mm_add_ps(lhs.m_fMat4x4[1], rhs.m_fMat4x4[1]);
	ret.m_fMat4x4[2] = _mm_add_ps(lhs.m_fMat4x4[2], rhs.m_fMat4x4[2]);
	ret.m_fMat4x4[3] = _mm_add_ps(lhs.m_fMat4x4[3], rhs.m_fMat4x4[3]);
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>         
gms::math::operator-(_In_ const SIMDMat4x4<gms::math::ArgsHolder<float>, float> &lhs, _In_ const SIMDMat4x4<gms::math::ArgsHolder<float>,float> &rhs){

	SIMDMat4x4<gms::math::ArgsHolder<float>,float> ret;
	ret.m_fMat4x4[0] = _mm_sub_ps(lhs.m_fMat4x4[0], rhs.m_fMat4x4[0]);
	ret.m_fMat4x4[1] = _mm_sub_ps(lhs.m_fMat4x4[1], rhs.m_fMat4x4[1]);
	ret.m_fMat4x4[2] = _mm_sub_ps(lhs.m_fMat4x4[2], rhs.m_fMat4x4[2]);
	ret.m_fMat4x4[3] = _mm_sub_ps(lhs.m_fMat4x4[3], rhs.m_fMat4x4[3]);
	return (ret);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>          
gms::math::operator*(_In_ const SIMDMat4x4<gms::math::ArgsHolder<float>, float> &lhs, _In_ SIMDMat4x4<gms::math::ArgsHolder<float>,float> &rhs){

	SIMDMat4x4<gms::math::ArgsHolder<float>,float> ret;
	rhs.transpose();
	__m128 t0(_mm_mul_ps(lhs.m_fMat4x4[0], rhs.m_fMat4x4[0]));
	ret.m_fMat4x4[0] = _mm_add_ps(t0, _mm_permute_ps(t0, 21));
	__m128 t1(_mm_mul_ps(lhs.m_fMat4x4[1], rhs.m_fMat4x4[1]));
	ret.m_fMat4x4[1] = _mm_add_ps(t1, _mm_permute_ps(t1, 21));
	__m128 t2(_mm_mul_ps(lhs.m_fMat4x4[2], rhs.m_fMat4x4[2]));
	ret.m_fMat4x4[2] = _mm_add_ps(t2, _mm_permute_ps(t2, 21));
	__m128 t3(_mm_mul_ps(lhs.m_fMat4x4[3], rhs.m_fMat4x4[3]));
	ret.m_fMat4x4[3] = _mm_add_ps(t3, _mm_permute_ps(t3, 21));
	return (ret);
}
	
gms::math::SIMDMat4x4<gms::math::ArgsHolder<float>, float>         
gms::math::operator/(_In_ const SIMDMat4x4<gms::math::ArgsHolder<float>,float> &rhs, _In_ const float s){

	SIMDMat4x4<gms::math::ArgsHolder<float>,float> ret;
	
	_ASSERT(s != 0.f);
		const float inv_s{ 1.f / s };
		ret.m_fMat4x4[0] = _mm_mul_ps(rhs.m_fMat4x4[0], _mm_set1_ps(inv_s));
		ret.m_fMat4x4[1] = _mm_mul_ps(rhs.m_fMat4x4[1], _mm_set1_ps(inv_s));
		ret.m_fMat4x4[2] = _mm_mul_ps(rhs.m_fMat4x4[2], _mm_set1_ps(inv_s));
		ret.m_fMat4x4[3] = _mm_mul_ps(rhs.m_fMat4x4[3], _mm_set1_ps(inv_s));
		return (ret);
}




gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>,double>::SIMDMat4x4(){

	this->m_dMat4x4[0] = _mm256_setzero_pd();
	this->m_dMat4x4[1] = _mm256_setzero_pd();
	this->m_dMat4x4[2] = _mm256_setzero_pd();
	this->m_dMat4x4[3] = _mm256_setzero_pd();
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ const ArgsHolder<double> &args) {

	this->m_dMat4x4[0] = _mm256_set_pd(args.a00, args.a01, args.a02, args.a03);
	this->m_dMat4x4[1] = _mm256_set_pd(args.a10, args.a11, args.a12, args.a13);
	this->m_dMat4x4[2] = _mm256_set_pd(args.a20, args.a21, args.a22, args.a23);
	this->m_dMat4x4[3] = _mm256_set_pd(args.a30, args.a31, args.a32, args.a33);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ const __m256d &v) {

	this->m_dMat4x4[0] = v;
	this->m_dMat4x4[1] = v;
	this->m_dMat4x4[2] = v;
	this->m_dMat4x4[3] = v;
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ const std::initializer_list<__m256d> &ilist) {

	this->m_dMat4x4[0] = ilist.begin()[0];
	this->m_dMat4x4[1] = ilist.begin()[1];
	this->m_dMat4x4[2] = ilist.begin()[2];
	this->m_dMat4x4[3] = ilist.begin()[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ const double s) {

	this->m_dMat4x4[0] = _mm256_set1_pd(s);
	this->m_dMat4x4[1] = _mm256_set1_pd(s);
	this->m_dMat4x4[2] = _mm256_set1_pd(s);
	this->m_dMat4x4[3] = _mm256_set1_pd(s);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ double(&ar)[16]){

	
	this->m_dMat4x4[0] = _mm256_set_pd(ar[0], ar[1], ar[2], ar[3]);
	this->m_dMat4x4[1] = _mm256_set_pd(ar[4], ar[5], ar[6], ar[7]);
	this->m_dMat4x4[2] = _mm256_set_pd(ar[8], ar[9], ar[10], ar[11]);
	this->m_dMat4x4[3] = _mm256_set_pd(ar[12], ar[13], ar[14], ar[15]);
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ float(&ar)[16]){

	constexpr int N{ 4 };
	for (int i{ 0 }; i < N; ++i){
		for (int j{ 0 }; j < N; ++j){
			this->m_dMat4x4[i].m256d_f64[j] = ar[i * N + j]; //non-efficient because of scalar access.
		}
	}
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ __m256d(&ar)[4]){

	this->m_dMat4x4[0] = ar[0];
	this->m_dMat4x4[1] = ar[1];
	this->m_dMat4x4[2] = ar[2];
	this->m_dMat4x4[3] = ar[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ const SIMDMat4x4 &rhs) {

	this->m_dMat4x4[0] = rhs.m_dMat4x4[0];
	this->m_dMat4x4[1] = rhs.m_dMat4x4[1];
	this->m_dMat4x4[2] = rhs.m_dMat4x4[2];
	this->m_dMat4x4[3] = rhs.m_dMat4x4[3];
}

gms::math::SIMDMat4x4<gms::math::ArgsHolder<double>, double>::SIMDMat4x4(_In_ SIMDMat4x4 &&rhs) {

	this->m_dMat4x4[0] = std::move(rhs.m_dMat4x4[0]);
	this->m_dMat4x4[1] = std::move(rhs.m_dMat4x4[1]);
	this->m_dMat4x4[2] = std::move(rhs.m_dMat4x4[2]);
	this->m_dMat4x4[3] = std::move(rhs.m_dMat4x4[3]);
}

