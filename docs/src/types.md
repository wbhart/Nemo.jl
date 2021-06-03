# Types in Nemo

Nemo is fully compatible with AbstractAlgebra.jl, but specialises implementations of
various commonly used rings with a highly optimised C implementation, provided by the
C libraries wrapped by Nemo.

Below, we give a list of all of the specialised types available in Nemo that implement
rings using a specialised C library. The types of elements of the respective rings and
other mathematical structures are given, and in parentheses we list the types of the
parent objects of the given rings and structures.

  - Flint
     - `fmpz` (`FlintIntegerRing`)
     - `fmpq` (`FlintRationalField`)
     - `nmod` (`NmodRing`)
     - `fq_nmod` (`FqNmodFiniteField`)
     - `fq` (`FqFiniteField`)
     - `padic` (`FlintPadicField`)
     - `qadic` (`FlintQadicField`)
     - `fmpz_poly` (`FmpzPolyRing`)
     - `fmpq_poly` (`FmpqPolyRing`)
     - `nmod_poly` (`NmodPolyRing`)
     - `fmpz_mod_poly` (`FmpzModPolyRing`)
     - `fq_poly` (`FqPolyRing`)
     - `fq_nmod_poly` (`FqNmodPolyRing`)
     - `fmpz_mpoly` (`FmpzMPolyRing`)
     - `fmpz_rel_series` (`FmpzRelSeriesRing`)
     - `fmpz_abs_series` (`FmpzAbsSeriesRing`)
     - `fmpq_rel_series` (`FmpqRelSeriesRing`)
     - `fmpq_abs_series` (`FmpqAbsSeriesRing`)
     - `fmpz_mod_rel_series` (`FmpzModRelSeriesRing`)
     - `fmpz_mod_abs_series` (`FmpzModAbsSeriesRing`)
     - `nmod_rel_series` (`NmodRelSeriesRing`)
     - `fq_nmod_rel_series` (`FqNmodRelSeriesRing`)
     - `fq_nmod_abs_series` (`FqNmodAbsSeriesRing`)
     - `fq_rel_series` (`FqRelSeriesRing`)
     - `fq_abs_series` (`FqAbsSeriesRing`)
     - `fmpz_mat` (`FmpzMatSpace`)
     - `fmpq_mat` (`FmpqMatSpace`)
     - `nmod_mat` (`NmodMatSpace`)
     - `fq_nmod_mat` (`FqNmodMatSpace`)
     - `fq_mat` (`FqMatSpace`)
     - `perm` (`SymmetricGroup`)

  - Antic
     - `nf_elem` (`AnticNumberField`)

  - Arb
     - `arb` (`ArbField`)
     - `acb` (`AcbField`)
     - `arb_poly` (`ArbPolyRing`)
     - `acb_poly` (`AcbPolyRing`)
     - `arb_mat` (`ArbMatSpace`)
     - `acb_mat` (`AcbMatSpace`)

  - Calcium

     - `qqbar` (`CalciumQQBarField`)
     - `ca` (`CalciumField`)
