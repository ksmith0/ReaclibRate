/** @file
 *  @author Karl Smith
 */

#include "ReaclibRate.hpp"

/** Constructor for charged particle reactions. Specifies the number of 
 *  resonances as well as the charge and reduced mass of the reactants.
 *
 *  For charged particle reactions the non-resonant set of parameters has a0 
 *  fixed based on S(0) (See ReaclibRate::SetSFactor). The a1 term is set based
 *  on the charge and mass of the reactants in the following form:
 *  \f[
 *    -4.2486 (Z_1^2 Z_2^2 \mu)^{1/3}
 *  \f]
 *  where \f$ Z_1 \f$ and \f$ Z_2 \f$ are the charge and \f$ \mu \f$ of the 
 *  reaction reactants.
 *  The a3 through a5 terms are allowed to float and the a6 term is set to -2/3.
 *
 *  For narrow resonances the a0 and a1 term are set based on the resonance 
 *  energy and strength (See ReaclibRate::SetResonance). The a2 through a5 
 *  terms are fixed to 0. Finally, a6 is set to -3/2.
 *
 *  \note Neutron induced non-resonant reaction rates are not yet supported.
 */
ReaclibRate::ReaclibRate(
	const char* name, const unsigned int numResonances, 
	const unsigned int z1, const unsigned int z2, const float mu
) :
	TF1(name, this, &ReaclibRate::Evaluate, 0.01, 10, 7 * (numResonances+1)), 
	numResonances_(numResonances),
	z1_(z1), 
	z2_(z2), 
	mu_amu_(mu) 
{
	//First set the non-resonant set of terms.
	//Set a0 as we do not yet know S(0).
	SetParameter(0, log(b_ * pow(z1_ * z2_ * mu_amu_, 1./3.)));
	FixParameter(1, 0);
	FixParameter(2, -4.2486 * pow(pow(z1_ * z2_, 2) * mu_amu_, 1./3.));
	//Parameters a3 through a5 are allowed to vary.
	FixParameter(6, -2./3.);

	//Now we set all resonant set terms.
	for (unsigned int i=0;i<numResonances_;i++) {
		//Set a0 as we do not yet know the strength.
		SetParameter(7 * (i+1) + 0, log(d_ * pow(mu_amu_, -3./2.)));
		//Set a1 as we do not yet know the resonance energy.
		SetParameter(7 * (i+1) + 1, -11.6045);
		for (int j=2;j<=5;j++) {
			FixParameter(7 * (i+1) + j, 0);
		}
		FixParameter(7 * (i+1) + 6, -3./2.);
	}
}

/** Sets the term (a0) of the non-resonant set associated with the s-factor at 
 *  energy zero, S(0). The a0 term takes the form 
 *  \f[
 *     ln[B (Z_1 Z_2 \mu)^{1/3} S(0)]
 *  \f] 
 *  where \f$ B = 7.8318 \times 10^9 cm^3 s^{-1} mole^{-1} MeV^{-1} \f$, 
 *  \f$ Z_1 \f$ and \f$ Z_2 \f$ are the charge of the reactants, \f$ \mu \f$ is
 *  the reduced mass of the reactants and \f$ S(0) \f$ is the value of the 
 *  S-factor evaluated at an energy of zero.
 *
 *  The parameter is fixed, but can be allowed to float by
 *  using the following: @code ReaclibRate::SetParLimits(0, 0, 0); @endcode
 */
void ReaclibRate::SetSFactor(float s0_MeVb) {
	FixParameter(0, log(b_ * pow(z1_ * z2_ * mu_amu_, 1./3.) * s0_MeVb));
}

/**Sets the terms for a resonance set. Specifically, a0 and a1 are set using 
 * the resonance strength and energy.
 *
 * For narrow resonances the a0 term takes on the following form:
 * \f[
 *   ln[D \mu^{-3/2} \omega\gamma]
 * \f]
 * where \f$ D = 1.5394 \times 10^{11} cm^3 s^{-1} mole^{-1} MeV^{-1} \f$,
 * \f$ \mu \f$ is the reactants reduced mass, and \f$ \omega \gamma \f$ is the
 * narrow resonance strength.
 *
 * The a1 term takes on the following form:
 * \f[
 *   -11.6045 E_r
 * \f]
 * where \f$ E_r \f$ is the resonance energy.
 *
 * The corresponding parameters are fixed, but can be allowed to float by
 * using the following: 
 * @code 
 * 	ReaclibRate::SetParLimits(7 * (resonanceId + 1) + 0, 0, 0);
 * 	ReaclibRate::SetParLimits(7 * (resonanceId + 1) + 1, 0, 0); 
 * @endcode
 */
void ReaclibRate::SetResonance(
	const unsigned int resonanceId, const float energy, const float strength
) {
	if (resonanceId <= numResonances_) {
		FixParameter(7 * (resonanceId+1) + 0, log(d_ * pow(mu_amu_, -3./2.) * strength));
		FixParameter(7 * (resonanceId+1) + 1, -11.6045 * energy);
	}
}

/**Extracts the reduced mass from the a2 term assuming that Z1 and Z2 are fixed.
 */
double ReaclibRate::GetReducedMass() {
	return pow(GetParameter(2) / -4.2486, 3.) / pow(z1_ * z2_, 2.);
}

/**Extracts the S-factor term at zero energy, S(0), from the term a0 using the 
 * reduced mass determined from a2.
 */
double ReaclibRate::GetSFactor() {
	float mu_amu = GetReducedMass();
	return exp(GetParameter(0)) / b_ / pow(z1_ * z2_ * mu_amu, 1./3.);
}

/**Extracts the resonance energy from the a1 term of the corresponding 
 * resonance set.
 */
double ReaclibRate::GetResonanceEnergy(const unsigned int resonanceId) {
	if (resonanceId > numResonances_) return -1;
	return GetParameter(7*resonanceId + 1) / -11.6045;
}

/**Extracts the resonance strength from the a0 term of the corresponding 
 * resonance set using the reduced mass determined from the a2 term of the 
 * non-resonant set.
 */
double ReaclibRate::GetResonanceStrength(const unsigned int resonanceId) {
	if (resonanceId > numResonances_) return -1;
	float mu_amu = GetReducedMass();
	return exp(GetParameter(7*resonanceId + 0)) / d_ / pow(mu_amu, -3./2.);
}

/**Evaluates the reaction by summing each set. The zeroth set is the 
 * non-resonant term while every additional set are resonant contributions.
 * The terms are evaluated using the following equation:
 * @f[
 * 	\sum_{n}exp\left[
 * 		a_{n,0} +\sum_{i=1}^5 a_{n,i} T_9^{2i-5/3} + a_{n,6} \ln T_9
 * 	\right]
 * @f]
 * where @f$ n @f$ is the set index. Typically a rate has one set for the 
 * non-resonant contribution and all subsequent sets are from narrow resonances.
 */
double ReaclibRate::Evaluate(double *t9, double *par) {
	//Value of the reaction rate at the specified t9 value.
	double reacRate = 0;
	//Loop over the resonances plus one for the non-resonant contribution.
	for (unsigned int i=0; i< numResonances_ + 1; i++) {
		//Compute the contribution from this component.
		// We apply the exponential at the end of each set.
		double component = par[7*i+0] + par[7*i+6] * log(t9[0]);
		for (int j=1;j<=5;j++) {
			component += par[7*i+j] * pow(t9[0],(2.*j-5.)/3.);
		}
		//Add the contribution to the total rate.
		reacRate += exp(component);
	}
	return reacRate;
}

