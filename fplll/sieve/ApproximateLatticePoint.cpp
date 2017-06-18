
template<class ET>
ApproximateLatticePoint<ET,-1>::ApproximateLatticePoint(ApproximateLatticePoint<ET,-1> const & other, int const dim)
    :length_exponent(other.length_exponent),approx(nullptr),approxn2(other.approxn2)
{
    approx = new ApproxEntryType[dim];
    memcpy(approx, other.approx, dim*sizeof(ApproxEntryType));
}

template <class ET> ApproximateLatticePoint<ET,-1>::ApproximateLatticePoint(ExactLatticePoint<ET,-1> const & exact_point)
{
    unsigned int const n = exact_point.get_dim(); //should be unsigned
    assert(n!=0);
    approx = new ApproxEntryType[n];
    static_assert (std::numeric_limits<ApproxEntryType>::is_specialized, "bad ApproxType");
    static_assert (std::numeric_limits<ApproxTypeNorm2>::is_specialized,"bad ApproxTypeNorm2");
    static_assert (std::numeric_limits<ApproxTypeNorm2>::digits /2 <= std::numeric_limits<ApproxEntryType>::digits, "bad Types"); //ensures that we have enough bits in ApproxEntryType. Note that /2 rounds down, which is correct.
    static_assert ( (! is_same<ET,Z_NR<mpz_t>>::value) || ( ( std::numeric_limits<ApproxTypeNorm2>::digits / 2 ) <= std::numeric_limits<long>::digits )," " ); //otherwise, mpz_get_si in the GMP class may not work correctly.
    unsigned int const max_bits = ( std::numeric_limits<ApproxTypeNorm2>::digits - floor(log2(n)) )/2;
    //may use __builtin_clz on GCC for improved efficiency, floor(log2(n)) is ridiculous...

    //max_bits now holds the number of bits we may use in the mantissas to prevent overflows when computing the norm2.

    signed int number_length = std::numeric_limits<signed int>::min(); //initial value, should be understood as -INFTY
    for(unsigned int i=0; i<n;++i)
    {
        number_length = std::max(number_length, LatticeApproximationsNew::get_exponent ( (exact_point)[i] ) );
    }
    //Now, number_length is the length of the longest coordinate (in absolute values);

    if(number_length == std::numeric_limits<signed int>::min()) // exact_point is all-zero vector. should never happen in the Sieve algorithm.
    {
        cerr << "Warning: approximating all-zero vector.";
        length_exponent=0;
        for(unsigned int i=0;i<n;++i){approx[i]=0;}
        approxn2=0;
        return;
    }
    else      // length_exponent is minimal, s.t. abs ( exact_point[i] / 2^length_exponent ) < 1 for all i
    {
        for(unsigned int i=0;i<n;++i)
        {
            length_exponent = number_length - max_bits;
            if(LatticeApproximations::MaybeRational<ET>::val == false) {length_exponent = max(0,length_exponent);} //constexpr if, actually...
            //The above line prevents padding with 0s from the right, if exact_point alreads fits into the mantissa.

            approx[i] = LatticeApproximationsNew::do_approximate<ApproxEntryType,ET> ( exact_point[i], length_exponent );
        }
        approxn2 = LatticeApproximationsNew::do_approximate<ApproxTypeNorm2,ET> ( exact_point.access_norm2(),2*length_exponent   );
    }
}
