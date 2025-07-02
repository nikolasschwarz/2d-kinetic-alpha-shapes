#include "Polynomial.hpp"

kinDS::Polynomial kinDS::operator+(const kinDS::CoefficientTerm& term1, const kinDS::CoefficientTerm& term2)
{
    return kinDS::Polynomial({ term1 }) + kinDS::Polynomial({ term2 });
}

kinDS::Polynomial kinDS::operator-(const kinDS::CoefficientTerm& term1, const kinDS::CoefficientTerm& term2)
{
    return kinDS::Polynomial({ term1 }) - kinDS::Polynomial({ term2 });
}
