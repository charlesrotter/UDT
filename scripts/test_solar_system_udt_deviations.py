#!/usr/bin/env python3
"""
Solar System UDT Deviations from GR
====================================

Tests for measurable deviations from General Relativity predictions
in the solar system due to finite R₀ effects in UDT.

This script calculates:
1. Orbital precession modifications
2. Light deflection differences
3. Gravitational redshift variations
4. Time delay modifications

These represent the most precise tests of UDT vs GR.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
import json
import os

class SolarSystemUDTTester:
    """Test UDT deviations from GR in solar system."""
    
    def __init__(self):
        """Initialize with solar system parameters."""
        # Physical constants
        self.c = 2.998e8           # Speed of light (m/s)
        self.G = 6.674e-11         # Gravitational constant
        self.AU = 1.496e11         # Astronomical unit (m)
        
        # Solar system masses
        self.M_sun = 1.989e30      # Solar mass (kg)
        self.M_earth = 5.972e24    # Earth mass (kg)
        self.M_jupiter = 1.898e27  # Jupiter mass (kg)
        
        # Characteristic scales
        self.rs_sun = 2 * self.G * self.M_sun / self.c**2  # Solar Schwarzschild radius
        
        # Planetary orbital data
        self.planets = {
            'Mercury': {'a': 0.387 * self.AU, 'e': 0.206, 'period': 87.97},
            'Venus': {'a': 0.723 * self.AU, 'e': 0.007, 'period': 224.7},
            'Earth': {'a': 1.000 * self.AU, 'e': 0.017, 'period': 365.25},
            'Mars': {'a': 1.524 * self.AU, 'e': 0.093, 'period': 687.0},
            'Jupiter': {'a': 5.203 * self.AU, 'e': 0.049, 'period': 4332.6}
        }
        
        # Results directory
        self.results_dir = "results/solar_system_udt_tests"
        os.makedirs(self.results_dir, exist_ok=True)
    
    def test_orbital_precession_deviations(self):
        """Test orbital precession differences between UDT and GR."""
        print("=" * 70)
        print("TESTING ORBITAL PRECESSION: UDT vs GR DEVIATIONS")
        print("=" * 70)
        print()
        
        # Range of UDT scales to test
        R0_values = [1e3 * self.AU, 1e4 * self.AU, 1e5 * self.AU, 1e6 * self.AU]
        
        precession_results = {}
        
        print("Planet    | Semi-major | GR Precession | UDT Precession (R₀=10³AU) | Deviation")
        print("          | Axis (AU)  | (\"/century)  | (\"/century)              | (\"/cent)")
        print("-" * 75)
        
        for planet, data in self.planets.items():
            a = data['a']  # Semi-major axis
            e = data['e']  # Eccentricity
            
            # GR precession formula
            precession_gr = (6 * np.pi * self.G * self.M_sun) / (self.c**2 * a * (1 - e**2))
            precession_gr_arcsec = precession_gr * (180/np.pi) * 3600 * 100  # arcsec/century
            
            # UDT precession for smallest R₀ (largest deviation)
            R0_test = R0_values[0]  # 1000 AU
            tau_perihelion = R0_test / (R0_test + a * (1 - e))  # Closest approach
            enhancement = 1 / tau_perihelion**2
            
            precession_udt = precession_gr * enhancement
            precession_udt_arcsec = precession_udt * (180/np.pi) * 3600 * 100
            
            deviation = precession_udt_arcsec - precession_gr_arcsec
            
            print(f"{planet:8s}  | {a/self.AU:8.3f}   | {precession_gr_arcsec:11.3f}   | {precession_udt_arcsec:19.3f}        | {deviation:+8.4f}")
            
            # Store results for all R₀ values
            planet_results = {
                'semi_major_axis_AU': a / self.AU,
                'eccentricity': e,
                'gr_precession_arcsec_per_century': precession_gr_arcsec,
                'udt_deviations': {}
            }
            
            for R0 in R0_values:
                tau_avg = R0 / (R0 + a)  # Average temporal dilation
                enhancement_avg = 1 / tau_avg**2
                precession_udt_r0 = precession_gr * enhancement_avg
                precession_udt_r0_arcsec = precession_udt_r0 * (180/np.pi) * 3600 * 100
                deviation_r0 = precession_udt_r0_arcsec - precession_gr_arcsec
                
                planet_results['udt_deviations'][f'R0_{R0/self.AU:.0f}_AU'] = {
                    'udt_precession_arcsec_per_century': precession_udt_r0_arcsec,
                    'deviation_from_gr_arcsec_per_century': deviation_r0,
                    'fractional_deviation': deviation_r0 / precession_gr_arcsec if precession_gr_arcsec > 0 else 0
                }
            
            precession_results[planet] = planet_results
        
        print()
        print("Observational Precision Requirements:")
        print("Mercury precession known to ~0.0001 arcsec/century")
        print("Venus precession measurements ~0.001 arcsec/century precision")
        print()
        
        return precession_results
    
    def test_light_deflection_modifications(self):
        """Test light deflection differences between UDT and GR."""
        print("=" * 70)
        print("TESTING LIGHT DEFLECTION: UDT vs GR MODIFICATIONS")
        print("=" * 70)
        print()
        
        # Impact parameters from solar limb to outer solar system
        impact_parameters = np.logspace(np.log10(self.rs_sun), np.log10(100*self.AU), 100)
        
        # UDT scales to test
        R0_values = [1e3 * self.AU, 1e4 * self.AU, 1e5 * self.AU]
        
        print("Testing light deflection for different R₀ values...")
        print()
        
        deflection_results = {}
        
        # GR light deflection: θ = 4GM/(c²b) for impact parameter b
        deflection_gr = 4 * self.G * self.M_sun / (self.c**2 * impact_parameters)
        deflection_gr_arcsec = deflection_gr * (180/np.pi) * 3600  # Convert to arcseconds
        
        print("Impact Parameter | GR Deflection | UDT Deflection (R₀=10³AU) | Deviation")
        print("(Solar Radii)    | (arcsec)      | (arcsec)                  | (μarcsec)")
        print("-" * 70)
        
        # Solar radius for reference
        R_sun = 6.96e8  # meters
        
        for i, b in enumerate(impact_parameters[::10]):  # Sample every 10th point
            # UDT modification
            R0_test = R0_values[0]  # 1000 AU
            tau = R0_test / (R0_test + b)
            enhancement = 1 / tau**2
            
            deflection_udt = deflection_gr[i*10] * enhancement
            deflection_udt_arcsec = deflection_udt * (180/np.pi) * 3600
            
            deviation_microarcsec = (deflection_udt_arcsec - deflection_gr_arcsec[i*10]) * 1e6
            
            print(f"{b/R_sun:12.1f}     | {deflection_gr_arcsec[i*10]:11.4f}   | {deflection_udt_arcsec:21.4f}    | {deviation_microarcsec:+9.2f}")
        
        # Store detailed results
        for R0 in R0_values:
            tau_array = R0 / (R0 + impact_parameters)
            enhancement_array = 1 / tau_array**2
            deflection_udt_array = deflection_gr * enhancement_array
            deflection_udt_arcsec_array = deflection_udt_array * (180/np.pi) * 3600
            
            deviation_array = deflection_udt_arcsec_array - deflection_gr_arcsec
            
            deflection_results[f'R0_{R0/self.AU:.0f}_AU'] = {
                'impact_parameters_m': impact_parameters.tolist(),
                'gr_deflection_arcsec': deflection_gr_arcsec.tolist(),
                'udt_deflection_arcsec': deflection_udt_arcsec_array.tolist(),
                'deviation_arcsec': deviation_array.tolist(),
                'max_deviation_microarcsec': np.max(np.abs(deviation_array)) * 1e6
            }
        
        print()
        print("Current observational precision: ~1 microarcsecond (Gaia, VLBI)")
        print()
        
        return deflection_results
    
    def test_gravitational_redshift_variations(self):
        """Test gravitational redshift modifications in UDT."""
        print("=" * 70)
        print("TESTING GRAVITATIONAL REDSHIFT: UDT vs GR VARIATIONS")
        print("=" * 70)
        print()
        
        # Distance range from solar surface to Earth orbit
        distances = np.logspace(np.log10(self.rs_sun), np.log10(self.AU), 100)
        
        # UDT scales
        R0_values = [1e3 * self.AU, 1e4 * self.AU, 1e5 * self.AU]
        
        redshift_results = {}
        
        print("Testing frequency shift for light escaping solar gravitational field...")
        print()
        
        # GR gravitational redshift: Δν/ν = -GM/(c²r)
        redshift_gr = -self.G * self.M_sun / (self.c**2 * distances)
        
        print("Distance     | GR Redshift   | UDT Redshift (R₀=10³AU) | Deviation")
        print("(Solar Radii)| (Δν/ν)        | (Δν/ν)                  | (parts in 10¹²)")
        print("-" * 70)
        
        R_sun = 6.96e8
        
        for i, r in enumerate(distances[::10]):
            R0_test = R0_values[0]
            tau = R0_test / (R0_test + r)
            
            # UDT redshift includes temporal dilation factor
            redshift_udt = redshift_gr[i*10] / tau**2
            
            deviation = (redshift_udt - redshift_gr[i*10]) * 1e12
            
            print(f"{r/R_sun:9.1f}    | {redshift_gr[i*10]:11.2e}   | {redshift_udt:21.2e}    | {deviation:+13.2f}")
        
        # Store detailed results
        for R0 in R0_values:
            tau_array = R0 / (R0 + distances)
            redshift_udt_array = redshift_gr / tau_array**2
            deviation_array = (redshift_udt_array - redshift_gr) * 1e12
            
            redshift_results[f'R0_{R0/self.AU:.0f}_AU'] = {
                'distances_m': distances.tolist(),
                'gr_redshift': redshift_gr.tolist(),
                'udt_redshift': redshift_udt_array.tolist(),
                'deviation_parts_per_trillion': deviation_array.tolist(),
                'max_deviation_ppt': np.max(np.abs(deviation_array))
            }
        
        print()
        print("Current precision: ~10⁻¹⁵ (atomic clocks, GPS)")
        print()
        
        return redshift_results
    
    def create_deviation_plots(self, precession_results, deflection_results, redshift_results):
        """Create comprehensive plots of UDT deviations from GR."""
        print("Creating solar system deviation plots...")
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Orbital precession deviations
        planets = list(precession_results.keys())
        R0_1000AU_deviations = [precession_results[p]['udt_deviations']['R0_1000_AU']['deviation_from_gr_arcsec_per_century'] 
                               for p in planets]
        orbital_distances = [precession_results[p]['semi_major_axis_AU'] for p in planets]
        
        ax1.semilogy(orbital_distances, np.abs(R0_1000AU_deviations), 'bo-', linewidth=2, markersize=8)
        for i, planet in enumerate(planets):
            ax1.annotate(planet, (orbital_distances[i], np.abs(R0_1000AU_deviations[i])), 
                        xytext=(5, 5), textcoords='offset points', fontsize=9)
        
        ax1.axhline(y=0.0001, color='r', linestyle='--', alpha=0.7, label='Current precision')
        ax1.set_xlabel('Orbital Distance (AU)')
        ax1.set_ylabel('|UDT - GR| Precession ("/century)')
        ax1.set_title('Orbital Precession Deviations (R₀=1000 AU)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Light deflection deviations
        R0_key = 'R0_1000_AU'
        impact_param = np.array(deflection_results[R0_key]['impact_parameters_m'])
        deviation = np.array(deflection_results[R0_key]['deviation_arcsec']) * 1e6  # microarcsec
        
        R_sun = 6.96e8
        ax2.loglog(impact_param/R_sun, np.abs(deviation), 'g-', linewidth=2)
        ax2.axhline(y=1, color='r', linestyle='--', alpha=0.7, label='Current precision (1 μas)')
        ax2.set_xlabel('Impact Parameter (Solar Radii)')
        ax2.set_ylabel('|UDT - GR| Light Deflection (μas)')
        ax2.set_title('Light Deflection Deviations (R₀=1000 AU)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Gravitational redshift deviations
        distances = np.array(redshift_results[R0_key]['distances_m'])
        redshift_dev = np.array(redshift_results[R0_key]['deviation_parts_per_trillion'])
        
        ax3.loglog(distances/R_sun, np.abs(redshift_dev), 'm-', linewidth=2)
        ax3.axhline(y=1e3, color='r', linestyle='--', alpha=0.7, label='Current precision (10⁻¹⁵)')
        ax3.set_xlabel('Distance from Sun (Solar Radii)')
        ax3.set_ylabel('|UDT - GR| Redshift (parts per 10¹²)')
        ax3.set_title('Gravitational Redshift Deviations (R₀=1000 AU)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Scale dependence of deviations
        R0_scales = [1e3, 1e4, 1e5]  # AU
        mercury_deviations = []
        deflection_deviations = []
        redshift_deviations = []
        
        for R0_val in R0_scales:
            # Mercury precession deviation
            merc_dev = precession_results['Mercury']['udt_deviations'][f'R0_{R0_val:.0f}_AU']['deviation_from_gr_arcsec_per_century']
            mercury_deviations.append(abs(merc_dev))
            
            # Light deflection at solar limb
            defl_key = f'R0_{R0_val:.0f}_AU'
            defl_dev = deflection_results[defl_key]['deviation_arcsec'][0] * 1e6  # First point, microarcsec
            deflection_deviations.append(abs(defl_dev))
            
            # Redshift at solar surface
            red_dev = redshift_results[defl_key]['deviation_parts_per_trillion'][0]
            redshift_deviations.append(abs(red_dev))
        
        ax4.loglog(R0_scales, mercury_deviations, 'bo-', label='Mercury precession ("/cent)', linewidth=2)
        ax4.loglog(R0_scales, deflection_deviations, 'gs-', label='Light deflection (μas)', linewidth=2)
        ax4.loglog(R0_scales, redshift_deviations, 'm^-', label='Redshift (ppt)', linewidth=2)
        
        ax4.set_xlabel('UDT Scale R₀ (AU)')
        ax4.set_ylabel('|UDT - GR| Deviation')
        ax4.set_title('Scale Dependence of UDT Deviations')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'{self.results_dir}/solar_system_udt_deviations.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Solar system deviation plots saved: {self.results_dir}/solar_system_udt_deviations.png")
        print()
    
    def assess_observational_detectability(self, results):
        """Assess which UDT deviations are potentially detectable."""
        print("=" * 70)
        print("OBSERVATIONAL DETECTABILITY ASSESSMENT")
        print("=" * 70)
        print()
        
        detectability = {}
        
        print("Effect                    | Best UDT Deviation | Current Precision | Detectable?")
        print("-" * 75)
        
        # Mercury precession (most precise)
        mercury_dev = abs(results['precession']['Mercury']['udt_deviations']['R0_1000_AU']['deviation_from_gr_arcsec_per_century'])
        precession_precision = 0.0001  # arcsec/century
        precession_detectable = mercury_dev > precession_precision
        
        print(f"Mercury precession        | {mercury_dev:.6f} \"/cent     | 0.0001 \"/cent    | {'YES' if precession_detectable else 'NO'}")
        
        # Light deflection at solar limb
        deflection_dev = results['deflection']['R0_1000_AU']['max_deviation_microarcsec']
        deflection_precision = 1.0  # microarcsec
        deflection_detectable = deflection_dev > deflection_precision
        
        print(f"Light deflection          | {deflection_dev:.2f} μas        | 1.0 μas           | {'YES' if deflection_detectable else 'NO'}")
        
        # Gravitational redshift
        redshift_dev = results['redshift']['R0_1000_AU']['max_deviation_ppt']
        redshift_precision = 1000  # parts per trillion (10^-15 = 1000 ppt)
        redshift_detectable = redshift_dev > redshift_precision
        
        print(f"Gravitational redshift    | {redshift_dev:.0f} ppt         | 1000 ppt          | {'YES' if redshift_detectable else 'NO'}")
        
        print()
        
        detectability = {
            'mercury_precession': {
                'deviation': mercury_dev,
                'precision': precession_precision,
                'detectable': precession_detectable,
                'confidence': 'HIGH' if mercury_dev > 5*precession_precision else 'LOW'
            },
            'light_deflection': {
                'deviation': deflection_dev,
                'precision': deflection_precision,
                'detectable': deflection_detectable,
                'confidence': 'HIGH' if deflection_dev > 5*deflection_precision else 'LOW'
            },
            'gravitational_redshift': {
                'deviation': redshift_dev,
                'precision': redshift_precision,
                'detectable': redshift_detectable,
                'confidence': 'HIGH' if redshift_dev > 5*redshift_precision else 'LOW'
            }
        }
        
        return detectability
    
    def run_solar_system_tests(self):
        """Run complete solar system UDT deviation tests."""
        print("\n" + "=" * 70)
        print("SOLAR SYSTEM UDT DEVIATION TESTS")
        print("=" * 70)
        print()
        
        print("Testing for measurable differences between UDT and GR...")
        print("in the solar system where R₀ is large but finite.")
        print()
        
        # Run tests
        precession_results = self.test_orbital_precession_deviations()
        deflection_results = self.test_light_deflection_modifications()
        redshift_results = self.test_gravitational_redshift_variations()
        
        # Combine results
        all_results = {
            'precession': precession_results,
            'deflection': deflection_results,
            'redshift': redshift_results
        }
        
        # Create plots
        self.create_deviation_plots(precession_results, deflection_results, redshift_results)
        
        # Assess detectability
        detectability = self.assess_observational_detectability(all_results)
        all_results['detectability'] = detectability
        
        # Save results
        with open(f'{self.results_dir}/solar_system_udt_deviations.json', 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        
        print("=" * 70)
        print("SOLAR SYSTEM TEST SUMMARY")
        print("=" * 70)
        print()
        print("UDT predicts measurable deviations from GR in solar system")
        print("if characteristic scale R₀ is in the range 10³-10⁴ AU.")
        print()
        print("MOST PROMISING TESTS:")
        print("1. Mercury orbital precession (highest precision)")
        print("2. Light deflection near solar limb (space-based measurements)")
        print("3. Gravitational redshift (atomic clock comparisons)")
        print()
        print("These tests could distinguish UDT from GR and")
        print("potentially measure the characteristic scale R₀!")
        print()
        print(f"Full results saved: {self.results_dir}/")
        
        return all_results

def main():
    """Main solar system testing function."""
    tester = SolarSystemUDTTester()
    results = tester.run_solar_system_tests()
    return results

if __name__ == "__main__":
    main()