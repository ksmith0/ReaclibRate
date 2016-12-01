/// @file
/// @author Karl Smith

#ifndef REACLIBRATE_H
#define REACLIBRATE_H

#include "TF1.h"

/**@brief A class inheriting from a TF1 that assists in the fitting of a 
 *   reaction rate to the JINA REACLIB format.
 * @author Karl Smith
 *
 * The REACLIB format is a set of seven parameters that provide a functional 
 * form for a reaction rate. The expression is designed to provide terms that
 * will fit both non-resonant, charge particle and neutron induced, as well as
 * narrow resonances. This class will take physical parameters (reactant 
 * charges, reduced mass, S(0), number of resonances, resonance energies, 
 * and resonance strengths) and make educated guess about initial values. After 
 * the fit is performed the resulting physical parameters can be extracted.
 */
class ReaclibRate : public TF1 {
	public:
		/// @brief Charged particle constructor.
		/// @param[in] name The name given to the rate.
		/// @param[in] numResonances The number of sets of resonance to add in 
		///   addition to the non-resonant set.
		/// @param[in] z1 The atomic number of the target.
		/// @param[in] z2 The atomic number of the reactant.
		/// @param[in] mu The reduced mass of the reactants in amu.
		ReaclibRate(
			const char* name, const unsigned int numResonances, 
			const unsigned int z1, const unsigned int z2, const float mu
		); 

		/// @brief Sets the best guess for the S-factor term S(0). 
		/// @param[in] s0_MeVb The value for the S-factor in MeV-b at energy zero, 
		///   S(0).
		void SetSFactor(const float s0_MeVb);

		/// @brief Set parameters for the specified narrow resonance.
		/// @param[in] resonanceId The ID of the associated resonance, begins at 0.
		/// @param[in] energy The resonance energy in MeV.
		/// @param[in] strength The resonance strength.
		void SetResonance(const unsigned int resonanceId, 
			const float energy, const float strength
		);

		/// @brief Evaluates the rate at the given temperature and parameters.
		/// @param[in] t9 Pointer to T9 values.
		/// @param[in] par Pointer to function parameters.
		/// @return The reaction rate for the specified t9 value and parameters.
		double Evaluate(double *t9, double *par);

		/// @brief Returns the S-factor, S(0), determined from the fit parameters.
		/// @return The S-factor in MeV-b.
		double GetSFactor();

		/// @brief Returns the reduced mass, determined from the fit parameters.
		/// @return The reduced mass in amu.
		double GetReducedMass();

		/// @brief Returns the resonance energy of the specified resonance set.
		/// @param[in] resonanceId The ID of the resonance to query, starts at 0.
		/// @return The resonance energy. If the resonance ID is invalid -1 is 
		///   returned.
		double GetResonanceEnergy(const unsigned int resosanceId);

		/// @brief Returns the resonance strength of the specified resonance set.
		/// @param[in] resonanceId The ID of the resonance to query, starts at 0.
		/// @return The resonance strength. If the resonance ID is invalid -1 is 
		///   returned.
		double GetResonanceStrength(const unsigned int resosanceId);

	private:
		const unsigned int numResonances_; ///< The number of resonance sets for this rate.
		const unsigned int z1_; ///< Atomic number of the target.
		const unsigned int z2_; ///< Atomic number of the reactant.
		const float mu_amu_; ///< Reduced mass of the reactants in amu.
		///Constant used for non-resonant a0 term. In units of @f$ cm^3 s^{-1} mole^{-1} MeV^{-1} barn^{-1} @f$
		static constexpr float b_ = 7.8318E9; 
		///Constant used for resonant a0 term. In units of @f$ cm^3 s^{-1} mole^{-1} MeV^{-1} @f$
		static constexpr float d_ = 1.5394E11;
};

#endif //REACLIBRATE_H
