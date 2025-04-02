#include <iostream>
#include <array>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <queue>
#include <mutex>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <immintrin.h>
#include <omp.h>
#include <csignal>
#include <random>
#include <algorithm>
#include <getopt.h>

#ifdef _WIN32
    #include <windows.h>
#else
    #include <unistd.h>
#endif

#include "sha256_avx2.h"
#include "ripemd160_avx2.h"
#include "SECP256K1.h"
#include "Point.h"
#include "Int.h"
#include "IntGroup.h"

using namespace std;

// Cross-platform terminal functions
void initConsole() {
#ifdef _WIN32
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    GetConsoleMode(hConsole, &mode);
    SetConsoleMode(hConsole, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
#endif
}

void clearTerminal() {
#ifdef _WIN32
    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    COORD coord = {0, 0};
    DWORD count;
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    GetConsoleScreenBufferInfo(hStdOut, &csbi);
    FillConsoleOutputCharacter(hStdOut, ' ', csbi.dwSize.X * csbi.dwSize.Y, coord, &count);
    SetConsoleCursorPosition(hStdOut, coord);
#else
    std::cout << "\033[2J\033[H";
#endif
    std::cout.flush();
}

void moveCursorTo(int x, int y) {
#ifdef _WIN32
    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    COORD coord = {(SHORT)x, (SHORT)y};
    SetConsoleCursorPosition(hStdOut, coord);
#else
    std::cout << "\033[" << y << ";" << x << "H";
#endif
    std::cout.flush();
}

// Configuration
int PUZZLE_NUM = 25;
int WORKERS = max(1, (int)thread::hardware_concurrency());
int FLIP_COUNT = -1;
const __uint128_t REPORT_INTERVAL = 100000;
static constexpr int POINTS_BATCH_SIZE = 512;
static constexpr int HASH_BATCH_SIZE = 8;

// Puzzle database
const unordered_map<int, tuple<int, string, string>> PUZZLE_DATA = {
        {20, {8, "b907c3a2a3b27789dfb509b730dd47703c272868",  "357535"}}, 
    {21, {9, "29a78213caa9eea824acf08022ab9dfc83414f56",  "863317"}},
    {22, {11, "7ff45303774ef7a52fffd8011981034b258cb86b", "1811764"}}, 
    {23, {12, "d0a79df189fe1ad5c306cc70497b358415da579e", "3007503"}},
    {24, {9, "0959e80121f36aea13b3bad361c15dac26189e2f",  "5598802"}},
    {25, {12, "2f396b29b27324300d0c59b17c3abc1835bd3dbb", "14428676"}},
    {26, {14, "bfebb73562d4541b32a02ba664d140b5a574792f", "33185509"}},
    {27, {13, "0c7aaf6caa7e5424b63d317f0f8f1f9fa40d5560", "54538862"}},
    {28, {16, "1306b9e4ff56513a476841bac7ba48d69516b1da", "111949941"}},
    {29, {18, "5a416cc9148f4a377b672c8ae5d3287adaafadec", "227634408"}},
    {30, {16, "d39c4704664e1deb76c9331e637564c257d68a08", "400708894"}},
    {31, {13, "d805f6f251f7479ebd853b3d0f4b9b2656d92f1d", "1033162084"}},
    {32, {14, "9e42601eeaedc244e15f17375adb0e2cd08efdc9", "2102388551"}},
    {33, {15, "4e15e5189752d1eaf444dfd6bff399feb0443977", "3093472814"}},
    {34, {16, "f6d67d7983bf70450f295c9cb828daab265f1bfa", "7137437912"}},
    {35, {19, "f6d8ce225ffbdecec170f8298c3fc28ae686df25", "14133072157"}},
    {36, {14, "74b1e012be1521e5d8d75e745a26ced845ea3d37", "20112871792"}},
    {37, {23, "28c30fb11ed1da72e7c4f89c0164756e8a021d",   "42387769980"}},
    {38, {21, "b190e2d40cfdeee2cee072954a2be89e7ba39364", "100251560595"}},
    {39, {23, "0b304f2a79a027270276533fe1ed4eff30910876", "146971536592"}},
    {40, {20, "95a156cd21b4a69de969eb6716864f4c8b82a82a", "323724968937"}},
    {41, {25, "d1562eb37357f9e6fc41cb2359f4d3eda4032329", "1003651412950"}},
    {42, {24, "8efb85f9c5b5db2d55973a04128dc7510075ae23", "1458252205147"}},
    {43, {19, "f92044c7924e5525c61207972c253c9fc9f086f7", "2895374552463"}},
    {44, {24, "80df54e1f612f2fc5bdc05c9d21a83aa8d20791e", "7409811047825"}},
    {45, {21, "f0225bfc68a6e17e87cd8b5e60ae3be18f120753", "15404761757071"}},
    {46, {24, "9a012260d01c5113df66c8a8438c9f7a1e3d5dac", "19996463086597"}},
    {47, {27, "f828005d41b0f4fed4c8dca3b06011072cfb07d4", "51408670348612"}},
    {48, {21, "8661cb56d9df0a61f01328b55af7e56a3fe7a2b2", "119666659114170"}},
    {49, {30, "0d2f533966c6578e1111978ca698f8add7fffdf3", "191206974700443"}},
    {50, {29, "de081b76f840e462fa2cdf360173dfaf4a976a47", "409118905032525"}},
    {51, {25, "ef6419cffd7fad7027994354eb8efae223c2dbe7", "611140496167764"}},
    {52, {27, "36af659edbe94453f6344e920d143f1778653ae7", "2058769515153876"}},
    {53, {26, "2f4870ef54fa4b048c1365d42594cc7d3d269551", "4216495639600700"}},
    {54, {30, "cb66763cf7fde659869ae7f06884d9a0f879a092", "6763683971478124"}},
    {55, {31, "db53d9bbd1f3a83b094eeca7dd970bd85b492fa2", "9974455244496707"}},
    {56, {31, "48214c5969ae9f43f75070cea1e2cb41d5bdcccd", "30045390491869460"}},
    {57, {33, "328660ef43f66abe2653fa178452a5dfc594c2a1", "44218742292676575"}},
    {58, {28, "8c2a6071f89c90c4dab5ab295d7729d1b54ea60f", "138245758910846492"}},
    {59, {30, "b14ed3146f5b2c9bde1703deae9ef33af8110210", "199976667976342049"}},
    {60, {31, "cdf8e5c7503a9d22642e3ecfc87817672787b9c5", "525070384258266191"}},
    {61, {25, "68133e19b2dfb9034edf9830a200cfdf38c90cbd", "1135041350219496382"}},
    {62, {35, "e26646db84b0602f32b34b5a62ca3cae1f91b779", "1425787542618654982"}},
    {63, {34, "ef58afb697b094423ce90721fbb19a359ef7c50e", "3908372542507822062"}},
    {64, {34, "3ee4133d991f52fdf6a25c9834e0745ac74248a4", "8993229949524469768"}},
    {65, {37, "52e763a7ddc1aa4fa811578c491c1bc7fd570137", "17799667357578236628"}},
    {66, {35, "20d45a6a762535700ce9e0b216e31994335db8a5", "30568377312064202855"}},
    {67, {31, "739437bb3dd6d1983e66629c5f08c70e52769371", "46346217550346335726"}},
    {68, {34, "e0b8a2baee1b77fc703455f39d51477451fc8cfc", "132656943602386256302"}} 
};

// Global state
vector<unsigned char> TARGET_HASH160_RAW(20);
string TARGET_HASH160;
Int BASE_KEY;
atomic<bool> stop_event(false);
mutex result_mutex;
queue<tuple<string, __uint128_t, int>> results;

// 128-bit counter with AVX
union AVXCounter {
    __m256i vec;
    uint64_t u64[4];
    __uint128_t u128[2];
    
    AVXCounter() : vec(_mm256_setzero_si256()) {}
    AVXCounter(__uint128_t value) { store(value); }
    
    void increment() {
        __m256i one = _mm256_set_epi64x(0, 0, 0, 1);
        vec = _mm256_add_epi64(vec, one);
        if (u64[0] == 0) {
            __m256i carry = _mm256_set_epi64x(0, 0, 1, 0);
            vec = _mm256_add_epi64(vec, carry);
        }
    }
    
    void add(__uint128_t value) {
        __m256i add_val = _mm256_set_epi64x(0, 0, value >> 64, value);
        vec = _mm256_add_epi64(vec, add_val);
        if (u64[0] < (value & 0xFFFFFFFFFFFFFFFFULL)) {
            __m256i carry = _mm256_set_epi64x(0, 0, 1, 0);
            vec = _mm256_add_epi64(vec, carry);
        }
    }
    
    __uint128_t load() const {
        return (static_cast<__uint128_t>(u64[1]) << 64) | u64[0];
    }
    
    void store(__uint128_t value) {
        u64[0] = static_cast<uint64_t>(value);
        u64[1] = static_cast<uint64_t>(value >> 64);
        u64[2] = 0;
        u64[3] = 0;
    }
};

AVXCounter total_checked_avx;
__uint128_t total_combinations = 0;
vector<string> g_threadPrivateKeys;
mutex progress_mutex;

// Performance tracking
atomic<uint64_t> globalComparedCount(0);
atomic<uint64_t> localComparedCount(0);
double globalElapsedTime = 0.0;
double mkeysPerSec = 0.0;
chrono::time_point<chrono::high_resolution_clock> tStart;

// Helper functions
string formatElapsedTime(double seconds) {
    int hrs = static_cast<int>(seconds) / 3600;
    int mins = (static_cast<int>(seconds) % 3600) / 60;
    int secs = static_cast<int>(seconds) % 60;
    ostringstream oss;
    oss << setw(2) << setfill('0') << hrs << ":"
        << setw(2) << setfill('0') << mins << ":"
        << setw(2) << setfill('0') << secs;
    return oss.str();
}

string to_string_128(__uint128_t value) {
    if (value == 0) return "0";
    char buffer[50];
    char* p = buffer + sizeof(buffer);
    *--p = '\0';
    while (value != 0) {
        *--p = "0123456789"[value % 10];
        value /= 10;
    }
    return p;
}

void signalHandler(int signum) {
    stop_event.store(true);
    cout << "\nInterrupt received, shutting down...\n";
}

class CombinationGenerator {
    int n, k;
    vector<int> current;
    
public:
    CombinationGenerator(int n, int k) : n(n), k(k), current(k) {
        if (k > n) k = n;
        for (int i = 0; i < k; ++i) current[i] = i;
    }

    static __uint128_t combinations_count(int n, int k) {
        if (k > n) return 0;
        if (k * 2 > n) k = n - k;
        if (k == 0) return 1;

        __uint128_t result = 1;
        for(int i = 1; i <= k; ++i) {
            result = result * (n - k + i) / i;
        }
        return result;
    }

    const vector<int>& get() const { return current; }
  
    bool next() {
        int i = k - 1;
        while (i >= 0 && current[i] == n - k + i) --i;
        if (i < 0) return false;
        ++current[i];
        for (int j = i + 1; j < k; ++j)
            current[j] = current[j-1] + 1;
        return true;
    }
  
    void unrank(__uint128_t rank) {
        __uint128_t total = combinations_count(n, k);
        if (rank >= total) {
            current.clear();
            return;
        }
        
        current.resize(k);
        __uint128_t remaining_rank = rank;
        int a = n;
        int b = k;
        __uint128_t x = (total - 1) - rank;
        
        for (int i = 0; i < k; i++) {
            a = largest_a_where_comb_a_b_le_x(a, b, x);
            current[i] = (n - 1) - a;
            x -= combinations_count(a, b);
            b--;
        }
    }

private:
    int largest_a_where_comb_a_b_le_x(int a, int b, __uint128_t x) const {
        while (a >= b && combinations_count(a, b) > x) a--;
        return a;
    }
};

void computeHash160Batch(uint8_t pubKeys[][33], uint8_t hashResults[][20], int count) {
    alignas(32) array<array<uint8_t, 64>, HASH_BATCH_SIZE> shaInputs;
    alignas(32) array<array<uint8_t, 32>, HASH_BATCH_SIZE> shaOutputs;
    alignas(32) array<array<uint8_t, 64>, HASH_BATCH_SIZE> ripemdInputs;

    auto byteswap_uint32 = [](uint32_t value) -> uint32_t {
    #if defined(_WIN32)
        return _byteswap_ulong(value);
    #else
        return __builtin_bswap32(value);
    #endif
    };

    for (int i = 0; i < count; i++) {
        memset(shaInputs[i].data(), 0, 64);
        memcpy(shaInputs[i].data(), pubKeys[i], 33);
        shaInputs[i][33] = 0x80;
        *reinterpret_cast<uint32_t*>(&shaInputs[i][60]) = byteswap_uint32(33 * 8);
    }

    uint8_t* inPtr[HASH_BATCH_SIZE];
    uint8_t* outPtr[HASH_BATCH_SIZE];
    for (int i = 0; i < HASH_BATCH_SIZE; i++) {
        inPtr[i] = shaInputs[i].data();
        outPtr[i] = shaOutputs[i].data();
    }

    sha256avx2_8B(inPtr[0], inPtr[1], inPtr[2], inPtr[3],
                  inPtr[4], inPtr[5], inPtr[6], inPtr[7],
                  outPtr[0], outPtr[1], outPtr[2], outPtr[3],
                  outPtr[4], outPtr[5], outPtr[6], outPtr[7]);

    for (int i = 0; i < count; i++) {
        memset(ripemdInputs[i].data(), 0, 64);
        memcpy(ripemdInputs[i].data(), shaOutputs[i].data(), 32);
        ripemdInputs[i][32] = 0x80;
        *reinterpret_cast<uint32_t*>(&ripemdInputs[i][60]) = byteswap_uint32(256);
    }

    for (int i = 0; i < HASH_BATCH_SIZE; i++) {
        inPtr[i] = ripemdInputs[i].data();
        outPtr[i] = hashResults[i].data();
    }

    ripemd160avx2::ripemd160avx2_32(
        inPtr[0], inPtr[1], inPtr[2], inPtr[3],
        inPtr[4], inPtr[5], inPtr[6], inPtr[7],
        outPtr[0], outPtr[1], outPtr[2], outPtr[3],
        outPtr[4], outPtr[5], outPtr[6], outPtr[7]
    );
}

void worker(Secp256K1* secp, int bit_length, int flip_count, int threadId, 
            __uint128_t start, __uint128_t end) {
    uint8_t localPubKeys[HASH_BATCH_SIZE][33];
    uint8_t localHashResults[HASH_BATCH_SIZE][20];
    int pointIndices[HASH_BATCH_SIZE];
    
    __m256i target16 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(TARGET_HASH160_RAW.data()));
    
    vector<Point> plusPoints(POINTS_BATCH_SIZE);
    vector<Point> minusPoints(POINTS_BATCH_SIZE);
    for (int i = 0; i < POINTS_BATCH_SIZE; i++) {
        Int tmp; tmp.SetInt32(i);
        plusPoints[i] = secp->ComputePublicKey(&tmp);
        minusPoints[i] = plusPoints[i];
        minusPoints[i].y.ModNeg();
    }

    vector<Int> deltaX(POINTS_BATCH_SIZE);
    IntGroup modGroup(POINTS_BATCH_SIZE);
    vector<Int> pointBatchX(2 * POINTS_BATCH_SIZE);
    vector<Int> pointBatchY(2 * POINTS_BATCH_SIZE);

    CombinationGenerator gen(bit_length, flip_count);
    gen.unrank(start);

    __uint128_t count = 0;
    while (!stop_event.load() && count < (end - start)) {
        Int currentKey;
        currentKey.Set(&BASE_KEY);
        
        const vector<int>& flips = gen.get();
        for (int pos : flips) {
            Int mask; mask.SetInt32(1); mask.ShiftL(pos);
            Int temp; temp.Set(&currentKey); temp.ShiftR(pos);
            temp.IsEven() ? currentKey.Add(&mask) : currentKey.Sub(&mask);
        }

        string keyStr = currentKey.GetBase16();
        keyStr = string(64 - keyStr.length(), '0') + keyStr;
        {
            lock_guard<mutex> lock(progress_mutex);
            g_threadPrivateKeys[threadId] = keyStr;
        }

        Point startPoint = secp->ComputePublicKey(&currentKey);
        Int startPointX = startPoint.x, startPointY = startPoint.y, startPointXNeg = startPointX;
        startPointXNeg.ModNeg();

        for (int i = 0; i < POINTS_BATCH_SIZE; i += 4) {
            deltaX[i].ModSub(&plusPoints[i].x, &startPointX);
            deltaX[i+1].ModSub(&plusPoints[i+1].x, &startPointX);
            deltaX[i+2].ModSub(&plusPoints[i+2].x, &startPointX);
            deltaX[i+3].ModSub(&plusPoints[i+3].x, &startPointX);
        }
        modGroup.Set(deltaX.data());
        modGroup.ModInv();

        for (int i = 0; i < POINTS_BATCH_SIZE; i += 4) {
            for (int j = 0; j < 4; j++) {
                // Plus points
                Int deltaY; deltaY.ModSub(&plusPoints[i+j].y, &startPointY);
                Int slope; slope.ModMulK1(&deltaY, &deltaX[i+j]);
                Int slopeSq; slopeSq.ModSquareK1(&slope);
                
                pointBatchX[i+j].Set(&startPointXNeg);
                pointBatchX[i+j].ModAdd(&slopeSq);
                pointBatchX[i+j].ModSub(&plusPoints[i+j].x);
                
                Int diffX; diffX.ModSub(&startPointX, &pointBatchX[i+j]);
                diffX.ModMulK1(&slope);
                
                pointBatchY[i+j].Set(&startPointY);
                pointBatchY[i+j].ModNeg();
                pointBatchY[i+j].ModAdd(&diffX);

                // Minus points
                deltaY.ModSub(&minusPoints[i+j].y, &startPointY);
                slope.ModMulK1(&deltaY, &deltaX[i+j]);
                slopeSq.ModSquareK1(&slope);
                
                pointBatchX[POINTS_BATCH_SIZE+i+j].Set(&startPointXNeg);
                pointBatchX[POINTS_BATCH_SIZE+i+j].ModAdd(&slopeSq);
                pointBatchX[POINTS_BATCH_SIZE+i+j].ModSub(&minusPoints[i+j].x);
                
                diffX.ModSub(&startPointX, &pointBatchX[POINTS_BATCH_SIZE+i+j]);
                diffX.ModMulK1(&slope);
                
                pointBatchY[POINTS_BATCH_SIZE+i+j].Set(&startPointY);
                pointBatchY[POINTS_BATCH_SIZE+i+j].ModNeg();
                pointBatchY[POINTS_BATCH_SIZE+i+j].ModAdd(&diffX);
            }
        }

        int localBatchCount = 0;
        for (int i = 0; i < 2*POINTS_BATCH_SIZE && localBatchCount < HASH_BATCH_SIZE; i++) {
            Point tempPoint;
            tempPoint.x.Set(&pointBatchX[i]);
            tempPoint.y.Set(&pointBatchY[i]);
            
            localPubKeys[localBatchCount][0] = tempPoint.y.IsEven() ? 0x02 : 0x03;
            for (int j = 0; j < 32; j++) {
                localPubKeys[localBatchCount][1 + j] = pointBatchX[i].GetByte(31 - j);
            }
            pointIndices[localBatchCount] = i;
            localBatchCount++;

            if (localBatchCount == HASH_BATCH_SIZE) {
                computeHash160Batch(localPubKeys, localHashResults, HASH_BATCH_SIZE);
                localComparedCount += HASH_BATCH_SIZE;

                for (int j = 0; j < HASH_BATCH_SIZE; j++) {
                    __m256i cand = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(localHashResults[j]));
                    int mask = _mm256_movemask_epi8(_mm256_cmpeq_epi8(cand, target16));
                    
                    if ((mask & 0x0F) == 0x0F) {
                        bool fullMatch = true;
                        for (int k = 0; k < 20; k++) {
                            if (localHashResults[j][k] != TARGET_HASH160_RAW[k]) {
                                fullMatch = false;
                                break;
                            }
                        }
                        
                        if (fullMatch) {
                            auto tEndTime = chrono::high_resolution_clock::now();
                            globalElapsedTime = chrono::duration<double>(tEndTime - tStart).count();
                            mkeysPerSec = (double)(globalComparedCount + localComparedCount) / globalElapsedTime / 1e6;
                            
                            Int foundKey; foundKey.Set(&currentKey);
                            int idx = pointIndices[j];
                            if (idx < POINTS_BATCH_SIZE) {
                                Int offset; offset.SetInt32(idx);
                                foundKey.Add(&offset);
                            } else {
                                Int offset; offset.SetInt32(idx - POINTS_BATCH_SIZE);
                                foundKey.Sub(&offset);
                            }
                            
                            string hexKey = foundKey.GetBase16();
                            hexKey = string(64 - hexKey.length(), '0') + hexKey;
                            
                            lock_guard<mutex> lock(result_mutex);
                            results.push(make_tuple(hexKey, total_checked_avx.load() + count, flip_count));
                            stop_event.store(true);
                            return;
                        }
                    }
                }
                
                localBatchCount = 0;
            }
        }

        if (!gen.next()) break;
        count++;

        if (count % 1000 == 0) {
            __uint128_t global_count = total_checked_avx.load() + count;
            if (global_count % REPORT_INTERVAL == 0) {
                lock_guard<mutex> lock(progress_mutex);
                double progress = min(100.0, (double)global_count / total_combinations * 100.0);
                moveCursorTo(0, 10 + threadId);
                cout << "Thread " << threadId << ": " << fixed << setprecision(2) << progress << "%";
                cout.flush();
            }
        }
    }

    total_checked_avx.add(count);
}

void printUsage(const char* programName) {
    cout << "Usage: " << programName << " [options]\n";
    cout << "Options:\n";
    cout << "  -p, --puzzle NUM    Puzzle number to solve (default: 20)\n";
    cout << "  -t, --threads NUM   Number of CPU cores to use (default: all)\n";
    cout << "  -f, --flips NUM     Override default flip count for puzzle\n";
    cout << "  -h, --help          Show this help message\n";
}

int main(int argc, char* argv[]) {
    signal(SIGINT, signalHandler);
    
    // Parse command line
    int opt;
    while ((opt = getopt(argc, argv, "p:t:f:h")) != -1) {
        switch (opt) {
            case 'p': PUZZLE_NUM = atoi(optarg); break;
            case 't': WORKERS = atoi(optarg); break;
            case 'f': FLIP_COUNT = atoi(optarg); break;
            case 'h': printUsage(argv[0]); return 0;
            default: printUsage(argv[0]); return 1;
        }
    }

    tStart = chrono::high_resolution_clock::now();

    // Initialize secp256k1 context
    Secp256K1 secp;
    secp.Init();
    
    // Load puzzle data
    auto puzzle_it = PUZZLE_DATA.find(PUZZLE_NUM);
    if (puzzle_it == PUZZLE_DATA.end()) {
        cerr << "Error: Invalid puzzle number\n";
        return 1;
    }
    
    auto [DEFAULT_FLIP_COUNT, TARGET_HASH160_HEX, PRIVATE_KEY_DECIMAL] = puzzle_it->second;
    if (FLIP_COUNT == -1) FLIP_COUNT = DEFAULT_FLIP_COUNT;
    TARGET_HASH160 = TARGET_HASH160_HEX;
    
    // Convert target hash
    for (int i = 0; i < 20; i++) {
        TARGET_HASH160_RAW[i] = stoul(TARGET_HASH160.substr(i * 2, 2), nullptr, 16);
    }
    
    // Set base key
    BASE_KEY.SetBase10(const_cast<char*>(PRIVATE_KEY_DECIMAL.c_str()));
    
    // Validate base key
    if (BASE_KEY.GetBitLength() > PUZZLE_NUM) {
        cerr << "Base key exceeds puzzle bit length!\n";
        return 1;
    }
    
    // Calculate total combinations
    total_combinations = CombinationGenerator::combinations_count(PUZZLE_NUM, FLIP_COUNT);
    
    // Display info
    string paddedKey = BASE_KEY.GetBase16();
    size_t firstNonZero = paddedKey.find_first_not_of('0');
    paddedKey = paddedKey.substr(firstNonZero);
    paddedKey = "0x" + paddedKey;

    cout << "=======================================\n";
    cout << "== Mutagen Puzzle Solver by Denevron ==\n";
    cout << "=======================================\n";    
    cout << "Puzzle: " << PUZZLE_NUM << " (" << PUZZLE_NUM << "-bit)\n";
    cout << "Target: " << TARGET_HASH160.substr(0, 10) << "..." << TARGET_HASH160.substr(30) << "\n";
    cout << "Base Key: " << paddedKey << "\n";
    cout << "Flip count: " << FLIP_COUNT << "\n";
    cout << "Total combinations: " << to_string_128(total_combinations) << "\n";
    cout << "Using " << WORKERS << " threads\n\n";

    // Validate combination count
    CombinationGenerator testGen(PUZZLE_NUM, FLIP_COUNT);
    __uint128_t actualCount = 0;
    while (testGen.next()) actualCount++;
    if (actualCount != total_combinations) {
        cerr << "ERROR: Combination count mismatch! Expected " 
             << to_string_128(total_combinations) << " got " 
             << to_string_128(actualCount) << "\n";
        return 1;
    }

    // Initialize threads
    g_threadPrivateKeys.resize(WORKERS, "0");
    vector<thread> threads;
    
    // Distribute work
    __uint128_t comb_per_thread = total_combinations / WORKERS;
    __uint128_t remainder = total_combinations % WORKERS;
    __uint128_t start = 0;
    
    for (int i = 0; i < WORKERS; i++) {
        __uint128_t end = start + comb_per_thread + (i < remainder ? 1 : 0);
        threads.emplace_back(worker, &secp, PUZZLE_NUM, FLIP_COUNT, i, start, end);
        start = end;
    }
    
    // Monitor progress
    while (!stop_event.load() && total_checked_avx.load() < total_combinations) {
        this_thread::sleep_for(chrono::milliseconds(100));
        
        globalElapsedTime = chrono::duration<double>(chrono::high_resolution_clock::now() - tStart).count();
        mkeysPerSec = (double)total_checked_avx.load() / globalElapsedTime / 1e6;
        
        double progress = min(100.0, (double)total_checked_avx.load() / total_combinations * 100.0);
        
        moveCursorTo(0, 8);
        cout << "Progress: " << fixed << setprecision(2) << progress << "%\n";
        cout << "Checked: " << to_string_128(total_checked_avx.load()) << " / " 
             << to_string_128(total_combinations) << "\n";
        cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";
        cout << "Elapsed: " << formatElapsedTime(globalElapsedTime) << "\n";
        
        for (int i = 0; i < WORKERS; i++) {
            moveCursorTo(0, 12 + i);
            cout << "Thread " << i << ": " << g_threadPrivateKeys[i];
        }
        cout.flush();
    }
    
    // Clean up
    for (auto& t : threads) {
        if (t.joinable()) t.join();
    }
    
    // Results
    globalElapsedTime = chrono::duration<double>(chrono::high_resolution_clock::now() - tStart).count();
    mkeysPerSec = (double)total_checked_avx.load() / globalElapsedTime / 1e6;
    
    if (!results.empty()) {
        auto [hex_key, checked, flips] = results.front();
        string compactHex = hex_key;
        firstNonZero = compactHex.find_first_not_of('0');
        compactHex = "0x" + compactHex.substr(firstNonZero);

        cout << "\n\nSOLUTION FOUND!\n";
        cout << "Private key: " << compactHex << "\n";
        cout << "Bit flips: " << flips << "\n";
        cout << "Checked " << to_string_128(checked) << " combinations\n";
        cout << "Time: " << fixed << setprecision(2) << globalElapsedTime << " seconds\n";
        cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";
        
        // Save solution
        ofstream out("puzzle_" + to_string(PUZZLE_NUM) + "_solution.txt");
        if (out) {
            out << hex_key;
            out.close();
            cout << "Saved to puzzle_" << PUZZLE_NUM << "_solution.txt\n";
        }
    } else {
        cout << "\n\nNo solution found. Checked all " 
             << to_string_128(total_checked_avx.load()) << " combinations\n";
        cout << "Time: " << fixed << setprecision(2) << globalElapsedTime << " seconds\n";
        cout << "Speed: " << fixed << setprecision(2) << mkeysPerSec << " Mkeys/s\n";
    }
    
    return 0;
}