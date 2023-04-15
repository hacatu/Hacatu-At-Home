#define _POSIX_C_SOURCE 202304L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdalign.h>
#include <math.h>

#include <crater/opts.h>
#include <crater/hash.h>
#include <primesieve.h>
#include <nut/factorization.h>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;

typedef struct [[gnu::packed]]{
	u64 magic;
	u64 shared_data_size;
	u128 n;
	u128 m;
	u64 a_ub;
	u64 a_dub;
	u64 p_ub;
	u64 num_primes;
	cr8r_hashtbl_t q_visited_schema;
	u64 primes[];
} spell_shared_hdr_t;

typedef struct{
	u128 n;
	u128 m;
	u64 a_ub;
	u64 a_dub;
	u64 p_ub;
	size_t num_primes;
	cr8r_hashtbl_t q_visited;
	const u64 *primes;
	int fd;
	u128 partial_result;
} spell_shared_t;

static inline u64 hash_u128(const cr8r_base_ft*, const void *_a){
	const u64 *a = (const u64*)_a;
	return a[0]^a[1];
}

static inline int cmp_u128(const cr8r_base_ft*, const void *_a, const void *_b){
	u128 a = *(const u128*)_a, b = *(const u128*)_b;
	if(a == b){
		return 0;
	}else if(a < b){
		return -1;
	}
	return 1;
}

static cr8r_hashtbl_ft hashft_u128 = {
	.base={.size=sizeof(u128)},
	.hash=hash_u128,
	.cmp=cmp_u128,
	.load_factor=0.5
};
static char buf[40];
static const spell_shared_hdr_t *shared_map;
static spell_shared_t shared_data;

static u64 calc_a_ub(i128 n){
	i128 x = cbrtl(n/2.l);
	while(1){
		i128 dx = (x*(x + 1)*(2*x + 1) - n) / (6*x*(x + 1) + 1);
		if(!dx){
			break;
		}
		x -= dx;
	}
	while(x*(x + 1)*(2*x + 1) <= n){
		++x;
	}
	return (u64)(x-1);
}

static u64 calc_a_dub(i128 a){
	i128 x = (sqrtl(a + 1) - 1)/2;
	while(4*x*(x + 1) <= a){
		++x;
	}
	return (u64)(x-1);
}

static bool init_shared_data(u128 n, u128 m){
	shared_data.n = n;
	shared_data.m = m;
	shared_data.a_ub = calc_a_ub(n);
	shared_data.a_dub = calc_a_dub(shared_data.a_ub);
	shared_data.p_ub = u64_nth_root(shared_data.a_ub + 1, 2);
	printf("init_shared_data(\n  n=%s, m=", cr8r_sprint_u128(buf, n));
	printf("%s,\n  a_ub=%"PRIu64", a_dub=%"PRIu64",\n p_ub=%"PRIu64"\n)\n\n", cr8r_sprint_u128(buf, m), shared_data.a_ub, shared_data.a_dub, shared_data.p_ub);
	shared_data.primes = primesieve_generate_primes(1, shared_data.p_ub + 1, &shared_data.num_primes, UINT64_PRIMES);
	if(!shared_data.primes || !cr8r_hash_init(&shared_data.q_visited, &hashft_u128, shared_data.p_ub)){
		fprintf(stderr, "\e[1;31mAllocation failed!\e[0m\n");
		free((void*)shared_data.primes);
		return false;
	}
	return true;
}

static bool save_shared_data(const char *filename){
	int fd = open(filename, O_CREAT | O_RDWR | O_EXCL, 0644);
	if(fd < 0){
		fprintf(stderr, "\e[1;31mCould not open file \"%s\": %s\e[0m\n", filename, strerror(errno));
		return false;
	}
	u128 dummy = 4;
	while(shared_data.q_visited.table_b){
		// look up a value not in the table over and over to force incremental resizing to finish if needed
		cr8r_hash_get(&shared_data.q_visited, &hashft_u128, &dummy);
	}
	u64 shared_data_size = offsetof(spell_shared_hdr_t, primes) + shared_data.num_primes*sizeof(u64);
	u64 qv_ta_offset = (shared_data_size + alignof(u128) - 1) & ~alignof(u128);
	u64 qv_ta_size = shared_data.q_visited.len_a*sizeof(u128);
	u64 qv_fa_offset = qv_ta_offset + qv_ta_size;
	u64 qv_fa_size = (((shared_data.q_visited.len_a - 1) >> 5) + 1)*sizeof(u64);
	shared_data_size = qv_fa_offset + qv_fa_size;
	if(ftruncate(fd, shared_data_size)){
		fprintf(stderr, "\e[1;31mCould not truncate (extend) file \"%s\": %s\e[0m\n", filename, strerror(errno));
		close(fd);
		return false;
	}
	void *pa = mmap(NULL, shared_data_size, PROT_WRITE, MAP_SHARED, fd, 0);
	if(pa == MAP_FAILED){
		fprintf(stderr, "\e[1;31mCould not mmap file \"%s\": %s\e[0m\n", filename, strerror(errno));
		close(fd);
		return false;
	}
	spell_shared_hdr_t *hdr = pa;
	hdr->magic = 0xDEADBEEF;
	hdr->shared_data_size = shared_data_size;
	hdr->n = shared_data.n;
	hdr->m = shared_data.m;
	hdr->a_ub = shared_data.a_ub;
	hdr->a_dub = shared_data.a_dub;
	hdr->p_ub = shared_data.p_ub;
	hdr->num_primes = shared_data.num_primes;
	hdr->q_visited_schema = (cr8r_hashtbl_t){
		.cap = shared_data.q_visited.cap,
		.len_a = shared_data.q_visited.len_a,
		.full = shared_data.q_visited.full
	};
	memcpy(hdr->primes, shared_data.primes, shared_data.num_primes*sizeof(u64));
	memcpy(pa + qv_ta_offset, shared_data.q_visited.table_a, qv_ta_size);
	memcpy(pa + qv_fa_offset, shared_data.q_visited.flags_a, qv_fa_size);
	shared_data.fd = fd;
	shared_map = hdr;
	return true;
}

static bool load_shared_data(const char *filename){
	int fd = open(filename, O_RDONLY);
	if(fd < 0){
		fprintf(stderr, "\e[1;31mCould not open file: %s\e[0m\n", strerror(errno));
		return false;
	}
	struct stat stat_buf = {};
	if(fstat(fd, &stat_buf)){
		fprintf(stderr, "\e[1;31mCould not stat file: %s\e[0m\n", strerror(errno));
		return false;
	}
	u64 file_size = stat_buf.st_size;
	void *pa = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
	if(pa == MAP_FAILED){
		fprintf(stderr, "\e[1;31mCould not mmap file: %s\e[0m\n", strerror(errno));
		close(fd);
		return false;
	}
	shared_map = pa;
	if(shared_map->magic != 0xDEADBEEF || shared_map->shared_data_size != file_size){
		fprintf(stderr, "\e[1;31mByte order / header size error!\e[0m\n");
		munmap(pa, file_size);
		close(fd);
		return false;
	}
	u64 shared_data_size = offsetof(spell_shared_hdr_t, primes) + shared_map->num_primes*sizeof(u64);
	u64 qv_ta_offset = (shared_data_size + alignof(u128) - 1) & ~alignof(u128);
	u64 qv_ta_size = shared_map->q_visited_schema.len_a*sizeof(u128);
	u64 qv_fa_offset = qv_ta_offset + qv_ta_size;
	shared_data.n = shared_map->n;
	shared_data.m = shared_map->m;
	shared_data.a_ub = shared_map->a_ub;
	shared_data.a_dub = shared_map->a_dub;
	shared_data.p_ub = shared_map->p_ub;
	shared_data.num_primes = shared_map->num_primes;
	shared_data.q_visited = (cr8r_hashtbl_t){
		.table_a = pa + qv_ta_offset,
		.flags_a = pa + qv_fa_offset,
		.cap = shared_map->q_visited_schema.cap,
		.len_a = shared_map->q_visited_schema.len_a,
		.full = shared_map->q_visited_schema.full
	};
	shared_data.primes = shared_map->primes;
	shared_data.fd = fd;
	printf("Loaded shared data:\n  n=%s, m=", cr8r_sprint_u128(buf, shared_data.n));
	printf("%s,\n  a_ub=%"PRIu64", a_dub=%"PRIu64", p_ub=%"PRIu64"\n", cr8r_sprint_u128(buf, shared_data.m), shared_data.a_ub, shared_data.a_dub, shared_data.p_ub);
	printf("  num_primes=%"PRIu64", q_count=%"PRIu64"\n", shared_data.num_primes, shared_data.q_visited.full);
	return true;
}

void free_resources(bool is_initial_process){
	if(munmap((void*)shared_map, shared_map->shared_data_size)){
		fprintf(stderr, "\e[1;31mmunmap failed: %s\e[0m\n", strerror(errno));
	}
	if(close(shared_data.fd)){
		fprintf(stderr, "\e[1;31mclose failed: %s\e[0m\n", strerror(errno));
	}
	if(is_initial_process){
		free((void*)shared_data.primes);
		cr8r_hash_destroy(&shared_data.q_visited, &hashft_u128);
	}
}

static inline void init_sqf(u64 *sqf, u64 *sqpr, u64 bucket_start, u64 bucket_end){
	for(u64 i = bucket_start; i < bucket_end; ++i){
		sqf[i-bucket_start] = i;
		sqpr[i-bucket_start] = 1;
	}
}

static inline bool s_pell_subrange(u64 a_start, u64 a_stop, bool is_initial_process){
	u128 bucket_size = 100*shared_data.p_ub;
	u64 *sqf = malloc(bucket_size*sizeof(u64));
	u64 *sqpr = malloc(bucket_size*sizeof(u64));
	//u64 max_ys = ceill(logl(sqrtl(8.l*n*n + 1) + sqrtl(2)*2*n)/logl(5 + 2*sqrtl(2)));
	u64 max_ys = 100;
	u128 *ys = malloc(max_ys*sizeof(u128));
	if(!sqf || !sqpr || !ys){
		fprintf(stderr, "\e[1;31mAllocation failed!\e[0m\n");
		free(sqf);
		free(sqpr);
		free(ys);
		return false;
	}
	printf("This subrange contains %.4f%% of the a values\n", (a_stop - a_start)*100./shared_data.a_ub);
	u128 tot4 = 0;
	//process the possible values of a in buckets.  since we need sqf(a*(a+1)),
	//we have the buckets overlap by 1 to avoid needing a special case for the first bucket
	u64 counter = 0;
	for(u64 bucket_start = a_start; bucket_start < a_stop + 1; bucket_start += bucket_size - 1){
		u64 bucket_end = bucket_start + bucket_size;
		if(bucket_end > a_stop + 1){
			bucket_end = a_stop + 1;
		}
		init_sqf(sqf, sqpr, bucket_start, bucket_end);
		//printf("bucket: [%"PRIu64", %"PRIu64")\n", bucket_start, bucket_end);
		for(size_t i = 0; i < shared_data.num_primes; ++i){
			u64 p = shared_data.primes[i];
			u64 p2 = p*p;
			for(u64 ppow = 1; !__builtin_mul_overflow(ppow, p2, &ppow) && ppow < bucket_end;){
				u64 m1 = bucket_start + ppow - 1;
				m1 -= m1%ppow;
				for(u64 idx = m1 - bucket_start; idx < bucket_end - bucket_start; idx += ppow){
					sqf[idx] /= p2;
					sqpr[idx] *= p;
				}
			}
		}
		u64 sqf_a = sqf[0];
		u64 sqpr_a = sqpr[0];
		for(u64 a = bucket_start; a < bucket_end - 1; ++a){
			if(++counter == 16000000){
				fprintf(stderr, "%.4f%%\n", (a - a_start)*100./(a_stop - a_start));
				counter = 0;
			}
			u128 tmp = sqf[a + 1 - bucket_start];
			u128 q = sqf_a*tmp;
			sqf_a = tmp;
			tmp = sqpr[a + 1 - bucket_start];
			u128 A = sqpr_a*tmp;
			sqpr_a = tmp;
			if(A*A*q != (u128)a*(a + 1)){
				fprintf(stderr, "\e[1;31mIncorrect data for %"PRIu64": A=%s", a, cr8r_sprint_u128(buf, A));
				fprintf(stderr, ", q=%s\e[0m\n", cr8r_sprint_u128(buf, q));
				return false;//this kills the process so it is ok to leak
			}
			if(cr8r_hash_get(&shared_data.q_visited, &hashft_u128, &q)){
				continue;
			}
			if(is_initial_process){
				int status = 0;
				cr8r_hash_insert(&shared_data.q_visited, &hashft_u128, &q, &status);
				if(status != 1){
					fprintf(stderr, "\e[1;31mq_visited is corrupted for %"PRIu64"!\e[0m\n", a);
					return false;//this kills the process so it is ok to leak
				}
			}
			u128 x1 = 2*a + 1;
			u128 y1 = A;
			u128 xk = x1*x1 + 4*q*y1*y1;
			u128 yk = 2*x1*y1;
			u128 M = 2*shared_data.n/(y1*q);
			u64 k = 1;
			ys[0] = y1;
			u128 term = 0;
			while(yk <= M){
				if(k == max_ys){
					fprintf(stderr, "\e[1;31mToo many ys for %"PRIu64"!\e[0m\n", a);
					return false;//this kills the process so it is ok to leak
				}
				ys[k++] = yk;
				tmp = x1*xk + 4*q*y1*yk;
				yk = x1*yk + y1*xk;
				xk = tmp;
			}
			//print(f"ys={ys}")
			//print(f"as={[str(y2a(y, q)) for y in ys]}")
			for(u64 i = 0; i < k - 1; ++i){
				u128 yi = ys[i];
				M = 2*shared_data.n/(yi*q);
				u128 yj = ys[i + 1];
				if(yj > M){
					break;
				}
				//print(f"  yi={yi}")
				u128 pterm = yj%shared_data.m;
				//print(f"    yj={yj}")
				for(u64 j = i + 2; j < k; ++j){
					yj = ys[j];
					if(yj > M){
						break;
					}
					//print(f"    yj={yj}")
					pterm = (pterm + yj)%shared_data.m;
				}
				term = (term + yi*pterm)%shared_data.m;
			}
			//printf("sum for q=%s: ", cr8r_sprint_u128(buf, q));
			//printf("%s\n", cr8r_sprint_u128(buf, term));
			tot4 = (tot4 + q*term)%shared_data.m;
		}
	}
	free(sqf);
	free(sqpr);
	free(ys);
	shared_data.partial_result = tot4*((shared_data.m + 1)/2)%shared_data.m;
	return true;
}

static bool s_pell_initial(u128 n, u128 modulus){
	if(!init_shared_data(n, modulus)){
		fprintf(stderr, "\e[1;31mAllocation of primes failed!\e[0m\n");
		return false;
	}
	printf("Processing initial subrange for a in [1, %"PRIu64")\n", shared_data.a_dub + 1);
	if(!s_pell_subrange(1, shared_data.a_dub + 1, true)){
		return false;
	}
	printf("q_count at exit: %"PRIu64"\n", shared_data.q_visited.full);
	return true;
}

int main(int argc, char **argv){
	u128 n, modulus;
	cr8r_opt options[] = {
		CR8R_OPT_HELP(options, "Compute S using Pell Equations\n"
			"Written by hacatu\n\n"),
		CR8R_OPT_GENERIC_OPTIONAL(&n, NULL, NULL, "upper bound for sqrt(T_a*T_b)"),
		CR8R_OPT_GENERIC_OPTIONAL(&modulus, NULL, NULL, "modulus for result"),
		CR8R_OPT_END()
	};
	bool parse_failed = false;
	bool is_initial_process = false;
	do{
		if(argc != 5){
			parse_failed = true;
			break;
		}else if(!strcmp("--initial", argv[1])){
			is_initial_process = true;
		}else if(!strcmp("--subrange", argv[1])){
			is_initial_process = false;
		}else{
			fprintf(stderr, "\e[1;31mUnknown operating mode\e[0m\n");
			parse_failed = true;
			break;
		}
	}while(0);
	for(u64 i = 3; i < 5; ++i){
		if(!cr8r_opt_parse_u128(options + i - 2, argv[i])){
			fprintf(stderr, "\e[1;31mCould not read argument %"PRIu64"\e[0m\n", i);
			parse_failed = true;
		}
	}
	if(parse_failed){
		fprintf(stderr, "\e[1;31mInvalid arguments, please call either like:\e[0m\n");
		fprintf(stderr, "\e[1;31m%s --initial <file>:str n:u128 modulus:u128\e[0m\n", argv[0]);
		fprintf(stderr, "\e[1;31  Which will initialize the shared data file \"file\" and\e[0m\n");
		fprintf(stderr, "\e[1;31  compute a partial sum for a in [1, a_dub]; OR\e[0m\n");
		fprintf(stderr, "\e[1;31m%s --subrange <file>:str a_start:u64 a_stop:u64\e[0m\n", argv[0]);
		fprintf(stderr, "\e[1;31m  Where a partial sum for a in [a_start, a_stop) will be performed.\e[0m\n");
		fprintf(stderr, "\e[1;31mIt is up to the user to combine the partial results at the end.\e[0m\n");
		exit(1);
	}
	if(is_initial_process){
		if(!s_pell_initial(n, modulus)){
			exit(1);
		}else if(!save_shared_data(argv[2])){
			exit(1);
		}
	}else{
		printf("Performing partial sum for a in [%"PRIu64", %"PRIu64")\n", (uint64_t)n, (uint64_t)modulus);
		if(!load_shared_data(argv[2])){
			exit(1);
		}else if(!s_pell_subrange(n, modulus, false)){
			exit(1);
		}
	}
	free_resources(is_initial_process);
	printf("partial_result: %s\n", cr8r_sprint_u128(buf, shared_data.partial_result));
}

