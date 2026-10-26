.PHONY: test2d fft_test_c2c fft_test_r2c timing halo_test io_test fft_test_r2c_mg

# Just build the examples
examples: test2d fft_test_c2c fft_test_r2c timing halo_test io_test fft_test_r2c_mg
	@echo "Built the examples"

test2d:
	cd test2d; $(MAKE) $@
fft_test_c2c:
	cd fft_test_c2c; $(MAKE) $@
fft_test_r2c:
	cd fft_test_r2c; $(MAKE) $@
timing:
	cd timing; $(MAKE) $@
halo_test:
	cd halo_test; $(MAKE) $@
io_test:
	cd io_test; $(MAKE) $@
fft_test_r2c_mg:
	cd fft_test_r2c_mg; $(MAKE) $@

# test all the examples (individual Makefiels should take care of updating)
basic_test:
	cd test2d; $(MAKE) $@
	cd fft_test_c2c; $(MAKE) $@
	cd fft_test_r2c; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@
	cd fft_test_r2c_mg; $(MAKE) $@

clean:
	cd test2d; $(MAKE) $@
	cd fft_test_c2c; $(MAKE) $@
	cd fft_test_r2c; $(MAKE) $@
	cd timing; $(MAKE) $@
	cd halo_test; $(MAKE) $@
	cd io_test; $(MAKE) $@
	cd fft_test_r2c_mg; $(MAKE) $@
