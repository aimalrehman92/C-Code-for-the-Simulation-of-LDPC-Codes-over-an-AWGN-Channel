# Code-for-the-Simulation-of-LDPC-Codes-on-AWGN-Channel
This repository contains C++ simulation of LDPC codes on AWGN channel with channel gains that follow Rayleigh distribution.

Brief discriptions of files:

**bpskrayleigh-main.cpp** --> This is the main file where we simulate BPSK transmission from source over an AWGN channel.
We read degree distribution of the optimized LDPC code in this file. The source message is encoded with LDPC code,
the graph of which is created in LDPCfnsBPSK.h file. The encoded ensemble goes through decoding process at the receiver
where Sum-Product-Algorithm is deployed. All the messages passed over the edges of the graph are in log-likelihood form.

**LDPCfnsBPSK.h** --> This header file contains all the functions required in the simulation of LDPC codes on the the AWGN
channel. From graph generation to the the function of Sum-Product-Algorithm to BER computation.

**filehandlingfns.h** --> This file contains functions for handling files

**.dat files** --> These files contain the degree distributions for the LDPC codes. They are read in the main file and
the corresponding decoding graphs are generated.
