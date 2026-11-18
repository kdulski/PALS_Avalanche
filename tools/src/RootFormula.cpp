#include "RootFormula.h"
#include "Emg.h"

#include <TFormula.h>
#include <TInterpreter.h>

#include <memory>
#include <stdexcept>

namespace pals {

FormulaFunction compileFormula(const std::string& expression, std::vector<std::string>& parameterNames)
{
    static bool declared = false;
    if (!declared) {
        if (!gInterpreter->Declare(palsEmgSource()))
            throw std::runtime_error("Cannot declare pals_emg in the ROOT interpreter");
        declared = true;
    }

    // shared: copies of the model (TF1 copies its functor) use the same compiled formula
    auto formula = std::make_shared<TFormula>("pals_model_formula", expression.c_str(), false);
    if (!formula->IsValid()) throw std::runtime_error("Cannot compile model.formula: " + expression);
    if (formula->GetNdim() > 1)
        throw std::runtime_error("model.formula may depend only on x (found " + std::to_string(formula->GetNdim()) +
                                 " variables)");

    parameterNames.clear();
    for (int i = 0; i < formula->GetNpar(); ++i) parameterNames.push_back(formula->GetParName(i));

    return [formula](double x, const double* params) { return formula->EvalPar(&x, params); };
}

} // namespace pals
