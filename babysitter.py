#!/usr/bin/env python3
from multiprocessing import Pool
import asyncio
import sys
import itertools as itt

"""
init_shared_data(
  n=100000000000000000000000000000000000, m=136101521,
  a_ub=368403149863, a_dub=303480,
 p_ub=606962
)
"""

a_start = 303481
a_stop = 368403149864
a_len = a_stop - a_start
num_subranges = 100

async def a_do_one_subrange(a_start, a_stop):
	proc = await asyncio.create_subprocess_shell(fr'time ./p833d --subrange 10e35.dat {a_start} {a_stop} 2>&1 | tee >(grep -vP "(^\d+\.\d+%$)|(^bucket: \[\d+, \d+\)$)" > {a_start}.{a_stop}.log)')
	await proc.wait()

def do_one_subrange(a_bounds):
	asyncio.run(a_do_one_subrange(*a_bounds))

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Expected two arguments:")
		print(f"{sys.argv[0]} i j")
		print("Where subranges [i, j), out of [0, 100) are to be processed")
		sys.exit(1)
	with Pool(2) as pool:
		for result in pool.imap_unordered(do_one_subrange, itt.pairwise(a_start + n*a_len//num_subranges for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1)), 1):
			pass

