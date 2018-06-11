#define __CL_ENABLE_EXCEPTIONS
#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<random>
#include<CL/cl.hpp>
#include<Windows.h>

#define D_INF 0x0FFFFFFF
#define USE_GPU

#ifdef USE_GPU
cl::CommandQueue *clqueue;
cl::Kernel *clkernel;
cl::Context *clcontext;

void gpumatmul(int *mata, int *matb, int *matc, int row, int col, int cr) {
	// create buffers on device (allocate space on GPU)
	cl::Buffer buffer_a(*clcontext, CL_MEM_READ_WRITE, sizeof(int) * row * cr);
	cl::Buffer buffer_b(*clcontext, CL_MEM_READ_WRITE, sizeof(int) * cr * col);
	cl::Buffer buffer_c(*clcontext, CL_MEM_READ_WRITE, sizeof(int) * row * col);

	// push write commands to queue
	clqueue->enqueueWriteBuffer(buffer_a, CL_TRUE, 0, sizeof(int) * row * cr, mata);
	clqueue->enqueueWriteBuffer(buffer_b, CL_TRUE, 0, sizeof(int) * cr * col, matb);

	clkernel->setArg(0, buffer_a);
	clkernel->setArg(1, buffer_b);
	clkernel->setArg(2, buffer_c);
	clkernel->setArg(3, row);
	clkernel->setArg(4, cr);
	clkernel->setArg(5, col);

	cl::Event qevent;
	clqueue->enqueueNDRangeKernel(*clkernel, cl::NullRange, cl::NDRange(row, col), cl::NDRange(1, 1), NULL, &qevent);
	qevent.wait();
	clqueue->enqueueReadBuffer(buffer_c, CL_TRUE, 0, sizeof(int) * row * col, matc);
	clqueue->flush();
}

#endif
class Matrix {
private:
	int row, col;
	int *matrix;
	static bool useGPU;
public:
	int id;
	Matrix() {
		matrix = NULL;
	}

	Matrix(int row, int col) {
		this->row = row;
		this->col = col;
		matrix = new int[row*col];
	}

	int& get(int x, int y) {
		return matrix[x*col + y];
	}

	friend std::ostream& operator<<(std::ostream& os, Matrix& mat) {
		for (int x = 0; x < mat.row; x++) {
			for (int y = 0; y < mat.col; y++) {
				os << mat.get(x,y) << "\t";
			}
			os << "\n";
		}
		return os;
	}

	~Matrix() {
		if (matrix != NULL) delete matrix;
		matrix = NULL;
	}

	//WARN: it modifies r-value(matrix)!!!
	Matrix& operator*(Matrix& mat) {
		if (col != mat.row) {
			std::cerr << "can not multiplication" << "\n";
			return mat;
		}
		//std::cout << id <<", " <<mat.id << std::endl;
		int *temp = new int[row*mat.col]();
#ifdef USE_GPU
		gpumatmul(this->matrix, mat.matrix, temp, row, mat.col, col);
#else
		for (int x = 0; x < row; x++) {
			for (int y = 0; y < mat.col; y++) {
				for (int i = 0; i < col; i++) {
					temp[x*mat.col + y] += get(x, i) * mat.get(i, y);
				}
			}
		}
#endif
		mat.row = row;
		if (mat.matrix != NULL)
			delete[] mat.matrix;
		mat.matrix = temp;
		return mat;
	}
	static void ramdomMatrix(Matrix& mat, int row, int col);
};


void Matrix::ramdomMatrix(Matrix& mat, int row, int col) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> udist(-9, 9);
	if (mat.matrix != NULL) delete[] mat.matrix;
	mat.matrix = new int[row*col];
	mat.row = row;
	mat.col = col;
	for (int i = 0; i < row*col; i++) {
		mat.matrix[i] = udist(gen);
	}

}

void ramdomd(int *d, int n) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> udist(3, 99);
	for (int i = 0; i < n; i++) {
		d[i] = udist(gen);
	}
}
int counter = 0;
int bfmin(const int *d, int *p, int i, int j, const int &n) {
	counter++;
	if (i == j) return 0;
	if (i == j + 1) return d[i] * d[i + 1] * d[i + 2];
	int min = bfmin(d, p, i, i, n) + bfmin(d, p, i + 1, j, n) + d[i] * d[i + 1] * d[j + 1];
	int vmin;
	p[i*n + j] = i;
	for (int x = i+1; x <= j - 1; x++) {
		vmin = bfmin(d, p, i, x, n) + bfmin(d, p, x + 1, j, n) + d[i] * d[x + 1] * d[j + 1];
		if (min > vmin) {
			min = vmin;
			p[i*n + j] = x;
		}
	}
	return min;
}

void printorder(const int *p, int i, int j, const int &n) {
	if (i == j) std::cout << "M" << i ;
	else {
		int k = p[i*n + j];
		std::cout << "(";
		printorder(p, i, k, n);
		printorder(p, k + 1, j, n);
		std::cout << ")";
	}
}

Matrix& mulorder(const int *p, int i, int j, const int &n, Matrix *mat) {
	if (i == j) return mat[i];
	else {
		int k = p[i*n + j];
		return mulorder(p, i, k, n, mat)*mulorder(p, k + 1, j, n, mat);
	}
}

void minmul(const int *d, int *m, int* p, const int &n) {
	int j, min[2];
	for (int diag = 1; diag < n; diag++) {
		for (int i = 0; i < n - diag; i++) {
			j = i + diag;
			min[0] = m[i*n + i] + m[(i + 1)*n + j] + d[i] * d[i + 1] * d[j + 1];
			m[i*n + j] = min[0];
			p[i*n + j] = i;
			for (int k = i + 1; k < j; k++) {
				min[1] = m[i*n + k] + m[(k + 1)*n + j] + d[i] * d[k + 1] * d[j + 1];
				if (min[0] > min[1]) {
					min[0] = min[1];
					m[i*n + j] = min[0];
					p[i*n + j] = k;
				}
			}
		}
	}
}

void isvalid(int *p, int *np, int n) {
	int *s1 = p, *s2 = np;
	for (int i = 0; i < n; i++) {
		if (*s1 != *s2){
			std::cout << "not good.\n";
			break;
		}
		s1++; s2++;
	}
	std::cout << "good.\n";
}

int main(int argc, char* argv[]) {
	cl_int err = CL_SUCCESS;
	LARGE_INTEGER StartingTime, EndingTime, ElapsedMicroseconds;
	LARGE_INTEGER Frequency;
	try {
#ifdef USE_GPU
		std::vector<cl::Platform> platforms;
		cl::Platform::get(&platforms);
		if (platforms.size() == 0) {
			std::cout << "Platform size 0\n";
			return -1;

		}
		cl::Platform default_platform = platforms[0];
		std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";
		cl_context_properties properties[] =
		{ CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[0])(), 0 };
		cl::Context context(CL_DEVICE_TYPE_GPU, properties);

		std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
		cl::Device default_device = devices[0];
		std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";

		std::string kernel_code =
			"   __kernel void mat_mul(global const int *A, global const int *B, global int *C, "
			"                          const int row, const int col, const int cr) {"
			"       int i, j;"
			""
			"       i = get_global_id(1);"
			"       j = get_global_id(0);"
			""
			"       for (int k=0; k<cr; k++)"
			"			C[i*col+j] += A[i*cr+k] * B[k*col+j];"
			"   }";
		/*
		std::string min_mul_kernel_code =
			"__kernel void min_mul(global const int *d, global int *m, global int *p, const int n) {"
			"	int diag = get_global_id(1) + 1;"
			"	if( diag >= n ) return;"
			"	int j, min[2];"
			"	int i = get_global_id(0);"
			"	if( i>=n-diag ) return;"
			"	j = i + diag;"
			"	min[0] = m[i*n + i] + m[(i + 1)*n + j] + d[i] * d[i + 1] * d[j + 1];"
			"	m[i*n + j] = min[0];"
			"	p[i*n + j] = i;"
			"	for (int k = i + 1; k < j; k++) {"
			"		min[1] = m[i*n + k] + m[(k + 1)*n + j] + d[i] * d[k + 1] * d[j + 1];"
			"		if (min[0] > min[1]) {"
			"			min[0] = min[1];"
			"			m[i*n + j] = min[0];"
			"			p[i*n + j] = k;"
			"		}"
			"	}"
			"}";*/
		
		std::string min_mul_kernel_code =
			"__kernel void min_mul(global const int *d, global int *m, global int *p, const int diag, const int n) {"
			"	int j, min[2];"
			"	int i = get_global_id(0);"
			"	if( i>=n-diag ) return;"
			"	j = i + diag;"
			"	min[0] = m[i*n + i] + m[(i + 1)*n + j] + d[i] * d[i + 1] * d[j + 1];"
			"	m[i*n + j] = min[0];"
			"	p[i*n + j] = i;"
			"	for (int k = i + 1; k < j; k++) {"
			"		min[1] = m[i*n + k] + m[(k + 1)*n + j] + d[i] * d[k + 1] * d[j + 1];"
			"		if (min[0] > min[1]) {"
			"			min[0] = min[1];"
			"			m[i*n + j] = min[0];"
			"			p[i*n + j] = k;"
			"		}"
			"	}"
			"}";

		cl::Program::Sources source(1,
			std::make_pair(kernel_code.c_str(), kernel_code.length()));
		cl::Program::Sources min_source(1,
			std::make_pair(min_mul_kernel_code.c_str(), min_mul_kernel_code.length()));
		cl::Program program_ = cl::Program(context, source);
		cl::Program min_program_ = cl::Program(context, min_source);
		program_.build(devices);
		min_program_.build(devices);
		cl::CommandQueue queue(context, default_device);
		cl::Kernel kernel(program_, "mat_mul");
		cl::Kernel min_kernel(min_program_, "min_mul");


		clcontext = &context;
		clkernel = &kernel;
		clqueue = &queue;
#endif
		bool safety = false;
		int n = -1;
		while (n < 2 || n > 16384){
			std::cout << "enter the count of matrix (2~16384): ";
			std::cin >> n;
		}
		if (n > 22) {
			std::cout << "n is big. brute force safety on.\n";
			safety = true;
		}
		int *w = new int[n*n];
		//int d[] = { 5, 2, 3, 4, 6, 7, 8, 3, 4, 5, 6, 2, 4, 5, 6, 6, 2, 4, 9, 5, 8, 6, 4 };
		//int d[] = { 200, 199, 199, 188, 88, 88, 85, 88, 87, 65, 54, 83, 72, 44, 32, 21, 13, 9, 6, 7, 8, 3, 4, 5, 6, 2 };
		int *d = new int[n + 1];
		ramdomd(d, n + 1);
		int *p = new int[n*n];
		int *np = new int[n*n];
		for (int i = 0; i<n*n; i++) {
			w[i] = 0;
			p[i] = -1;
			np[i] = -1;
		}
		ULONGLONG delta;
		if (!safety) {
			std::cout << "brute force:\n";
			QueryPerformanceFrequency(&Frequency);
			QueryPerformanceCounter(&StartingTime);
			bfmin(d, np, 0, n-1, n);
			QueryPerformanceCounter(&EndingTime);
			//std::cout << counter << std::endl;
			ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
			ElapsedMicroseconds.QuadPart *= 1000000;
			ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;

			printf("ended in %llu.%06llu\n", ElapsedMicroseconds.QuadPart / 1000000, ElapsedMicroseconds.QuadPart);
		}
#ifdef USE_GPU
		cl::Buffer buffer_d(context, CL_MEM_READ_WRITE, sizeof(int) * (n + 1));
		cl::Buffer buffer_m(context, CL_MEM_READ_WRITE, sizeof(int) * n * n);
		cl::Buffer buffer_p(context, CL_MEM_READ_WRITE, sizeof(int) * n * n);
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartingTime);
		queue.enqueueWriteBuffer(buffer_d, CL_TRUE, 0, sizeof(int) * (n + 1), d);
		queue.enqueueWriteBuffer(buffer_m, CL_TRUE, 0, sizeof(int) * n * n, w);
		queue.enqueueWriteBuffer(buffer_p, CL_TRUE, 0, sizeof(int) * n * n, p);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		printf("mem cpy ended in %llu.%06llu\n", ElapsedMicroseconds.QuadPart / 1000000, ElapsedMicroseconds.QuadPart);

		std::cout << "dynamic programming(GPU):\n";
		
		min_kernel.setArg(0, buffer_d);
		min_kernel.setArg(1, buffer_m);
		min_kernel.setArg(2, buffer_p);
		min_kernel.setArg(4, n);
		int worksize = 256;
		while (n > worksize) {
			worksize += 256;
		}
		//delta = GetTickCount64();
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartingTime);
		for(int i=1;i<n;i++){
			min_kernel.setArg(3, i);
			queue.enqueueNDRangeKernel(min_kernel, cl::NullRange, cl::NDRange(worksize), cl::NDRange(16), NULL);
		}

		queue.enqueueReadBuffer(buffer_p, CL_TRUE, 0, sizeof(int) * n * n, p);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
		//delta = GetTickCount64() - delta;
		
#else
		std::cout << "dynamic programming:\n";
		QueryPerformanceFrequency(&Frequency);
		QueryPerformanceCounter(&StartingTime);
		minmul(d, w, np, n);
		QueryPerformanceCounter(&EndingTime);
		ElapsedMicroseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
		ElapsedMicroseconds.QuadPart *= 1000000;
		ElapsedMicroseconds.QuadPart /= Frequency.QuadPart;
#endif
		/*
		isvalid(p, np, n);
		std::cout << std::endl;
		*/
		//printf("ended in %llu.%03llu\n", delta / 1000, delta % 1000);
		printf("ended in %llu.%06llu\n", ElapsedMicroseconds.QuadPart / 1000000, ElapsedMicroseconds.QuadPart);
		std::cout << "generating random matrix..:\n";
		Matrix *mat = new Matrix[n]();
		for(int i=0;i<n;i++){
			mat[i].id = i;
			Matrix::ramdomMatrix(mat[i], d[i], d[i+1]);
			//std::cout << mat[i];

		}
		if (!safety) {
			std::cout << "*** brtue force optimal checking ***\n";
			std::cout << "brute force result:\n";
			printorder(np, 0, n - 1, n);
			std::cout << std::endl;
			std::cout << "dynamic programming result:\n";
			printorder(p, 0, n-1, n);
			std::cout << std::endl;
		}
		std::cout << "optimal multiplication:\n";
		delta = GetTickCount64();
		Matrix &result = mulorder(p, 0, n-1, n, mat);
		delta = GetTickCount64() - delta;
		printf("ended in %llu.%03llu\n", delta / 1000, delta % 1000);
		//std::cout << result;
		
		for (int i = 0; i<n; i++) {
			Matrix::ramdomMatrix(mat[i], d[i], d[i + 1]);
			//std::cout << mat[i];
		}
		std::cout << "non-optimal multiplication:\n";
		delta = GetTickCount64();
		for (int i = 0; i < n-1; i++) {
			mat[i] * mat[i + 1];
		}
		delta = GetTickCount64() - delta;
		printf("ended in %llu.%03llu\n", delta / 1000, delta % 1000);
		
		//std::cout << mat[0]*mat[1]*mat[2];
		//printorder(p, 0, 5);
		/*
		std::cout << std::endl;
		Matrix tet1(5, 5);
		Matrix tet2(5, 5);
		for (int x = 0; x < 5; x++) {
			for (int y = 0; y < 5; y++) {
				tet1.get(x, y) = x+y;
				tet2.get(x, y) = x+y;
			}
		}
		std::cout << tet1 * tet2 << std::endl;
		*/
		//printf("ended in %d.%03d", delta / 1000, delta % 1000);

	}
	catch (cl::Error err) {
		std::cerr
			<< "ERROR: "
			<< err.what()
			<< "("
			<< err.err()
			<< ")"
			<< std::endl;

	}
	std::cout << "press enter to exit...";
	getchar();
	getchar();
	return EXIT_SUCCESS;
}