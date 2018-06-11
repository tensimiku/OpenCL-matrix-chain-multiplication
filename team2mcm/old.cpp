/*
std::cout << "result: {\n";
for (int i = 0; i < n*n; i++) {
std::cout << w[i] << " ";
if (i % n == n - 1) std::cout << "\n";
}
for (int i = 0; i<n*n; i++) {
std::cout << p[i] << " ";
if (i % n == n - 1) std::cout << "\n";
}
for (int i = 0; i<n*n; i++) {
std::cout << np[i] << " ";
if (i % n == n - 1) std::cout << "\n";
}
*/
#ifdef USE_OLD
struct leaf {
	int head, tail;
	bool headsmall;
	int sup_weight;
	int cost;
	int base_product, side_product;
	struct leaf *desendant, *next;
};
int cp[MSIZE];
struct leaf *list;
int n;

void init_hushing(const int *d) {
	cp[0] = 0;
	for (int i = 1; i < MSIZE; i++) {
		cp[i] = cp[i - 1] + d[i - 1] * d[i];
	}
	list = new struct leaf[MSIZE + 1]();
}
void hs_stack_push(int &top_element, int &second_element, int &tos, int *stack, int c) {
	stack[tos - 1] = c;
	second_element = top_element;
	top_element = c;
	tos -= 1;
};
void hs_stack_pop(int &top_element, int &second_element, int &tos, int *stack) {
	tos += 1;
	top_element = second_element;
	second_element = stack[tos - 3];
};

void hushing_onesweep(const int *d) {

	int current, tos = MSIZE;
	int top_element = -1, second_element = -1;
	struct leaf *arc_list = list;
	int stack[MSIZE];
	hs_stack_push(top_element, second_element, tos, stack, 0);
	hs_stack_push(top_element, second_element, tos, stack, 1);
	current = 2;

	while (current < MSIZE) {
		if ((d[second_element] <= d[top_element]) && (d[top_element] > d[current])) {
			arc_list[current].head = second_element;
			arc_list[current].tail = current;
			arc_list[current].headsmall = d[second_element] <= d[current];
			arc_list[current].base_product = d[second_element] * d[current];
			arc_list[current].side_product = cp[current] - d[second_element];
			arc_list[current].desendant = NULL;
			arc_list[current].next = &arc_list[current + 1];

			hs_stack_pop(top_element, second_element, tos, stack);
			if (tos >= MSIZE - 1) {
				hs_stack_push(top_element, second_element, tos, stack, current);
				current += 1;
			}
		}
		else {
			hs_stack_push(top_element, second_element, tos, stack, current);
			current += 1;
		}
		while ((tos <= MSIZE - 3) && (d[second_element] <= d[top_element]) && (d[top_element] > d[MSIZE - 1])) {
			arc_list[current].head = second_element;
			arc_list[current].tail = MSIZE - 1;
			arc_list[current].headsmall = d[second_element] <= d[MSIZE - 1];
			arc_list[current].base_product = d[second_element] * d[MSIZE - 1];
			arc_list[current].side_product = cp[MSIZE - 1] - cp[second_element];
			arc_list[current].desendant = NULL;
			arc_list[current].next = &arc_list[current + 1];
		}
	}
}
#endif
/*
int bf_hu_shing(const int *d, int i, int j, const std::vector<int> &visible, int n) {
if (j == i + 2)
return d[i] * d[i + 1] * d[j];
if (j<i + 2)
return 0;
int min = i+1;
int to;
std::vector<int> s_visible; // for recursive

for (int k = i + 1; k <= j; k++) {
if (d[min] > d[k]) {
min = k;
}
}
to = (min + 1) % n;
for (int k = i+1; k <= j-1; k++) {
if (min == k) continue;
if (d[to] > d[k]) {
to = k;
}
}
if (min > to) {
for (int k = i; k <= j; k++) {
if (k <= to || k >= min) {
s_visible.push_back(k);
}
}
int half = reversed_hu_shing(d, min, to, s_visible, s_visible.size());
s_visible.clear();
for (int k = i; k <= j; k++) {
if (k >= to || k <= min) {
s_visible.push_back(k);
}
}
half += bf_hu_shing(d, to, min, s_visible, n - s_visible.size());
return half;
}
for (int k = i; k <= j; k++) {
if (k >= to || k <= min) {
s_visible.push_back(k);
}
}
int half = bf_hu_shing(d, min, to, s_visible, n - s_visible.size());
s_visible.clear();
for (int k = i; k <= j; k++) {
if (k <= to || k >= min) {
s_visible.push_back(k);
}
}
half += reversed_hu_shing(d, to, min, s_visible, s_visible.size());
return half;
}

int hu_shing(const int *d, int i, int j, const std::vector<int> &visible, int n) {
std::cout << "n value :" << n << std::endl;
if (n == 3) {
std::cout << "returned from i, j:" << i << ", " << j << std::endl;
return d[visible.at(0)] * d[visible.at(1)] * d[visible.at(2)];
}
if (n < 3)
return 0;
int min = j;
int to;
std::vector<int> s_visible; // for recursive

for (int k = 0; k < n; k++) {
to = visible.at(k);
if (to != min && to != i )
break;
}
for (int k = 0; k < n; k++) {
int idx = visible.at(k);
if (i == idx || min == idx) continue;
if (d[to] > d[idx]) {
to = idx;
}
}
if (min > to) {
for (int k = 0; k < n; k++) {
int idx = visible.at(k);
if (idx <= to || idx >= min) {
s_visible.push_back(idx);
}
}
int half = hu_shing(d, min, to, s_visible, s_visible.size());
s_visible.clear();
for (int k = i; k <= j; k++) {
int idx = visible.at(k);
if (idx >= to || idx <= min) {
s_visible.push_back(idx);
}
}
half += hu_shing(d, to, min, s_visible, s_visible.size());
return half;
}
for (int k = 0; k < n; k++) {
int idx = visible.at(k);
if (idx >= to || idx <= min) {
s_visible.push_back(k);
}
}
int half = hu_shing(d, min, to, s_visible, s_visible.size());
s_visible.clear();
for (int k = i; k <= j; k++) {
if (k <= to || k >= min) {
s_visible.push_back(k);
}
}
half += hu_shing(d, to, min, s_visible, s_visible.size());
return half;

}
int sweepstart(const int *d, int i, int j) {
if (j == i + 2)
return d[i] * d[i + 1] * d[j];
if (j<i + 2)
return 0;
int min = i;
int to;
std::vector<int> visible; // for reversed

for (int k = i + 1; k <= j; k++) {
if (d[min] > d[k]) {
min = k;
}
}
to = (min + 1)%(j-i+1);
for (int k = i; k <= j; k++) {
if (min == k) continue;
if (d[to] > d[k]) {
to = k;
}
}
for (int k = i; k <= j; k++) {
if (k <= to || k >= min) {
visible.push_back(k);
}
}
int half = hu_shing(d, min, to, visible, visible.size());
visible.clear();
for (int k = i; k <= j; k++) {
if (k >= to || k <= min) {
visible.push_back(k);
}
}
half += hu_shing(d, to, min, visible, visible.size());
return half;
}
*/