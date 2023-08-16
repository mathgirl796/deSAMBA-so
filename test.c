#include <stdio.h>
#include "src/lib/khash.h"

KHASH_MAP_INIT_INT(shit, void*)

int main() {
	int ret, is_missing;
	khiter_t k;
	khash_t(shit) *h = kh_init(shit);
	k = kh_put(shit, h, 88, &ret);
    printf("%d\t%p\n", kh_key(h, k), kh_val(h, k));
	kh_value(h, k) = (void*)666;
	k = kh_get(shit, h, 88);
    printf("%d\t%p\n", kh_key(h, k), kh_val(h, k));
	is_missing = (k == kh_end(h));
	k = kh_get(shit, h, 5);
    printf("%d\t%p\n", kh_key(h, k), kh_val(h, k));
	kh_del(shit, h, k);
	for (k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h, k)) kh_value(h, k) = (void*)1;
	kh_destroy(shit, h);
	return 0;
}