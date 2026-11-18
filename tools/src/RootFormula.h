#ifndef PALS_ROOTFORMULA_H
#define PALS_ROOTFORMULA_H

#include "FormulaModel.h"

#include <string>
#include <vector>

namespace pals {

// Compiles the formula from the config with TFormula (variable x, named
// parameters [name]). Besides the usual TFormula / TMath functions the formula
// may use pals_emg(x, sigma, t0, tau): exponential convolved with a Gaussian,
// unit area. Returns the function and the parameter names in TFormula order.
FormulaFunction compileFormula(const std::string& expression, std::vector<std::string>& parameterNames);

} // namespace pals

#endif
