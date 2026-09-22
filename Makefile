# Top-level targets.  All commands run in the foreground.
.PHONY: pdf verify-cubic verify-third smoke-third clean

pdf:
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

# Full regeneration of every cubic certificate log (~3 h on one core; set JOBS=k to parallelise the Laplace grid).
verify-cubic:
	cd cubic && ./run_all.sh

# Full regeneration of every Third-Borwein certificate log (several hours; see third/run_all.sh for SKIP_EXACT/SKIP_F).
verify-third:
	cd third && ./run_all.sh

# Fast subset (under two minutes): saddle endpoint, saddle coverage, piece (N), Laplace tail cells, band checks.
smoke-third:
	cd third/certificates && python3 centre3_saddle_end.py && python3 cls012_cov.py > /dev/null && echo "cls012_cov OK" && \
	python3 offpeak_tb.py 821 0.25 | tail -1 && python3 cls012_lap.py tail 11.99 | cut -c1-120 && \
	python3 lap34_cert.py tail 11.99 | cut -c1-120 && python3 bandfix_coeffs.py | tail -1

clean:
	cd paper && latexmk -c
