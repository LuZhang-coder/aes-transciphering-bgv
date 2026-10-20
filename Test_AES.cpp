#include <NTL/tools.h>
#include <NTL/BasicThreadPool.h>
/** Test_AES.cpp - adapted for HElib 2.2.0 */
#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <thread>
#include "homAES.h"
helib::SecKey* gVerifyKey = nullptr;  // 临时验证用（全局作用域）

using namespace std;
using namespace NTL;

static long mValues[][14] = {
  { 2,  512,    771, 16,771,  0,  0,     5,    0,    0,-32,  0,  0, 100},
  { 2, 4096,   4369, 16, 17, 257, 0,   258, 4115,    0, 16,-16,  0, 100},
  { 2, 16384, 21845, 16, 17, 5, 257,  8996,17477,21591, 16,  4,-16,1600},
  { 2, 23040, 28679, 24, 17, 7, 241, 15184, 4098,28204, 16,  6,-10,1500},
  { 2, 46080, 53261, 24, 17,13, 241, 43863,28680,15913, 16, 12,-10, 100},
  { 2, 64512, 65281, 48, 97,673,  0, 43073,22214,    0, 96,-14,  0, 100}
};

void printState(Vec<uint8_t>& st);
extern long AESKeyExpansion(unsigned char RoundKey[], unsigned char Key[], int NN);
extern void Cipher(unsigned char out[16], unsigned char in[16],
                   unsigned char RoundKey[], int Nr);

// 辅助函数：打印vector<long>
static void printVec(const vector<long>& v) {
  cout << "[";
  for (size_t i=0; i<v.size(); i++) { if(i) cout<<" "; cout<<v[i]; }
  cout << "]";
}
static void printVec(const Vec<long>& v) {
  cout << "[";
  for (long i=0; i<v.length(); i++) { if(i) cout<<" "; cout<<v[i]; }
  cout << "]";
}

extern int getInverseValue(int num);
extern int getSBoxValue(int num);
int main(int argc, char **argv)
{
  NTL::SetNumThreads(64);
  std::cout << "NTL threads set to: " << NTL::AvailableThreads() << std::endl;
  long idx = 0;
  long c   = 3;
  long L   = 0;
  long B   = 23;
  bool boot = false;

  // 简单命令行解析：key=value 格式
  for (int i=1; i<argc; i++) {
    string s(argv[i]);
    if      (s.rfind("sz=",0)==0)   idx  = stol(s.substr(3));
    else if (s.rfind("c=",0)==0)    c    = stol(s.substr(2));
    else if (s.rfind("L=",0)==0)    L    = stol(s.substr(2));
    else if (s.rfind("B=",0)==0)    B    = stol(s.substr(2));
    else if (s.rfind("boot=",0)==0) boot = (s.substr(5)=="1"||s.substr(5)=="true");
  }
  if (idx > 5) idx = 5;

  if (boot) {
    if (L < 23) L = 23;
    if (idx < 1) idx = 1;
  } else {
#if (NTL_SP_NBITS<50)
    if (L < 10) L = 10;
#else
    if (L < 10) L = 10;
#endif
  }

  long p = mValues[idx][0];
  long m = mValues[idx][2];

  Vec<long> mvec;
  vector<long> gens, ords;

  append(mvec, mValues[idx][4]);
  if (mValues[idx][5] > 1) append(mvec, mValues[idx][5]);
  if (mValues[idx][6] > 1) append(mvec, mValues[idx][6]);

  gens.push_back(mValues[idx][7]);
  if (mValues[idx][8] > 1) gens.push_back(mValues[idx][8]);
  if (mValues[idx][9] > 1) gens.push_back(mValues[idx][9]);

  ords.push_back(mValues[idx][10]);
  if (abs(mValues[idx][11]) > 1) ords.push_back(mValues[idx][11]);
  if (abs(mValues[idx][12]) > 1) ords.push_back(mValues[idx][12]);

  cout << "*** Test_AES: c=" << c << ", L=" << L << ", B=" << B
       << ", boot=" << boot << ", m=" << m << " (=";
  printVec(mvec);
  cout << "), gens="; printVec(gens);
  cout << ", ords="; printVec(ords);
  cout << endl;

  helib::setTimersOn();
  double tm = -NTL::GetTime();
  cout << "computing key-independent tables..." << flush;

  helib::Context context =
    helib::ContextBuilder<helib::BGV>()
      .m(m).p(p).r(1)
      .gens(gens).ords(ords)
      .bits(L * B)
      .c(c)
      .bootstrappable(boot)
      .mvec(mvec)
      .thickboot()
      .build();

  tm += NTL::GetTime();
  cout << "done in " << tm << " seconds\n";

  cout << "  " << context.numPrimes() << " primes, total bitsize=" << context.bitSizeOfQ()
       << ", security=" << context.securityLevel() << endl;
  cout << "  " << context.getZMStar().getNSlots() << " slots ("
       << (context.getZMStar().getNSlots()/16) << " blocks) per ctxt\n";

  cout << "computing key-dependent tables..." << flush;
  tm = -NTL::GetTime();

  helib::SecKey secretKey(context);
  helib::PubKey& publicKey = secretKey;
  secretKey.GenSecKey();
  gVerifyKey = &secretKey;  // 临时验证用

  long ord = context.getZMStar().OrderOf(0);
  for (long i = 1; i < 16; i++) {
    long exp = i * ord / 16;
    long val = NTL::PowerMod(context.getZMStar().ZmStarGen(0), exp, m);
    secretKey.GenKeySWmatrix(1, val);
    if (!context.getZMStar().SameOrd(0))
      secretKey.GenKeySWmatrix(1, NTL::InvMod(val, m));
  }
  helib::addFrbMatrices(secretKey);
  if (boot) {
    helib::addSome1DMatrices(secretKey);
    secretKey.genRecryptData();
  }

  tm += NTL::GetTime();
  cout << "done in " << tm << " seconds\n";

  cout << "computing AES tables..." << flush;
  tm = -NTL::GetTime();
  HomAES hAES(context);
  const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2 = hAES.getEA();
  cout << "[诊断] ea2.getDegree() = " << ea2.getDegree() << endl;
  long blocksPerCtxt = ea2.size() / 16;
  long nBlocks = blocksPerCtxt;

  Vec<uint8_t> ptxt(INIT_SIZE, nBlocks*16);
  Vec<uint8_t> aesCtxt(INIT_SIZE, nBlocks*16);
  Vec<uint8_t> aesKey(INIT_SIZE, 16);
  uint8_t keySchedule[240];

  { GF2X rnd;
    random(rnd, 8*ptxt.length());
    BytesFromGF2X(ptxt.data(), rnd, ptxt.length());
    random(rnd, 8*aesKey.length());
    BytesFromGF2X(aesKey.data(), rnd, aesKey.length()); }

  vector<helib::Ctxt> encryptedAESkey;
  hAES.encryptAESkey(encryptedAESkey, aesKey, publicKey);
  vector<vector<helib::Ctxt>> encryptedAESkeyA;
  hAES.encryptAESkey_singleBlock16Ctxt(encryptedAESkeyA, aesKey, publicKey);
  tm += NTL::GetTime();
  cout << "done in " << tm << " seconds\n";






  // ===== 方案A 已移除（已验证与方案B结果一致，保留方案B作为唯一实验）=====
  {
    // ===== [方案B] 1920-block 批量SIMD验证测试 =====
    cout << "\n\n######## [MARKER] >>>>> 方案B 开始 <<<<< ########\n" << endl;
    cerr << "\n\n######## [MARKER] >>>>> 方案B 开始 <<<<< ########\n" << endl;
    cout << "\n=== [方案B] " << ea2.size() << "-block 批量(SIMD batching) 测试 ===" << endl;
    {
      long numBlocks = ea2.size();
      Vec<uint8_t> multiBlocks(INIT_SIZE, numBlocks*16);
      { GF2X rnd;
        random(rnd, 8*multiBlocks.length());
        BytesFromGF2X(multiBlocks.data(), rnd, multiBlocks.length()); }

      vector<ZZX> multiEncoded;
      vector<vector<helib::Ctxt>> encryptedAESkeyB;
      hAES.encryptAESkey_broadcastBlock16Ctxt(encryptedAESkeyB, aesKey, publicKey);

      encode4AES_multiBlock16Ctxt(multiEncoded, multiBlocks.data(), numBlocks, ea2);
      cout << "[诊断] multiEncoded.size() = " << multiEncoded.size()
           << ", numBlocks = " << numBlocks << endl;

      vector<helib::Ctxt> eDataB;
      eDataB.reserve(16);
      for (long j=0; j<16; j++) {
        eDataB.push_back(helib::Ctxt(publicKey));
        publicKey.Encrypt(eDataB[j], multiEncoded[j]);
      }

      double tmB = -NTL::GetTime();
      hAES.bitLevelRoundTestA(eDataB, encryptedAESkeyB);
      tmB += NTL::GetTime();
      cout << "[方案B] 10轮总耗时(" << numBlocks << "个block, 不含最初加密) = "
           << tmB << " 秒 (平均每block " << (tmB/numBlocks) << " 秒)" << endl;

      {
        int correctCount = 0, incorrectCount = 0;
        for (long i=0; i<16; i++) {
          if (eDataB[i].isCorrect()) correctCount++;
          else incorrectCount++;
        }
        cout << "[FHE-CORRECTNESS][方案B] isCorrect() 检验: "
             << correctCount << "/16 通过, " << incorrectCount << "/16 未通过" << endl;
      }

      vector<ZZX> polyB(16);
      for (long j=0; j<16; j++) secretKey.Decrypt(polyB[j], eDataB[j]);
      Vec<uint8_t> multiCtxtBlocks(INIT_SIZE, numBlocks*16);
      decode4AES_multiBlock16Ctxt(multiCtxtBlocks.data(), numBlocks, polyB, ea2);

      unsigned char keyScheduleB[240];
      AESKeyExpansion(keyScheduleB, aesKey.data(), 128);

      long blockMatchCount = 0, blockMismatchCount = 0;
      for (long i=0; i<numBlocks; i++) {
        uint8_t refBlockB[16];
        Cipher(refBlockB, &multiBlocks[i*16], keyScheduleB, 10);
        bool blockMatch = true;
        for (int b=0; b<16; b++) {
          if (multiCtxtBlocks[i*16+b] != refBlockB[b]) { blockMatch = false; break; }
        }
        if (blockMatch) blockMatchCount++;
        else {
          blockMismatchCount++;
          if (blockMismatchCount <= 5)
            cerr << "[方案B] block " << i << " 不匹配!" << endl;
        }
      }
      cout << "[方案B] 逐block正确性: " << blockMatchCount << "/" << numBlocks
           << " 匹配, " << blockMismatchCount << "/" << numBlocks << " 不匹配" << endl;
      if (blockMismatchCount == 0)
        cout << "[方案B] 结果: 全部一致 ✓" << endl;
      else
        cerr << "[方案B] 结果: 存在不一致 ✗" << endl;
    }
    // ===== [方案B] 测试结束 =====
    cout << "\n\n######## [MARKER] >>>>> 方案B 结束 <<<<< ########\n" << endl;
    cerr << "\n\n######## [MARKER] >>>>> 方案B 结束 <<<<< ########\n" << endl;

    // ===== [方案C] 4组x16密文=64密文, 7680-block 批量测试 =====
    cout << "\n\n######## [MARKER] >>>>> 方案C 开始 <<<<< ########\n" << endl;
    cerr << "\n\n######## [MARKER] >>>>> 方案C 开始 <<<<< ########\n" << endl;
    {
      long numGroups = 4;
      long numBlocksPerGroup = ea2.size();
      long totalBlocksC = numGroups * numBlocksPerGroup;
      cout << "\n=== [方案C] " << numGroups << "组 x " << numBlocksPerGroup
           << "-block = " << totalBlocksC << "-block 批量(64密文) 测试 ===" << endl;

      Vec<uint8_t> multiBlocksC(INIT_SIZE, totalBlocksC*16);
      { GF2X rnd;
        random(rnd, 8*multiBlocksC.length());
        BytesFromGF2X(multiBlocksC.data(), rnd, multiBlocksC.length()); }

      vector<ZZX> multiEncodedC;
      vector<vector<helib::Ctxt>> encryptedAESkeyC;
      hAES.encryptAESkey_multiGroupBroadcastCtxt(encryptedAESkeyC, aesKey, numGroups, publicKey);

      encode4AES_multiGroupBlock16Ctxt(multiEncodedC, multiBlocksC.data(), numGroups, numBlocksPerGroup, ea2);
      cout << "[诊断] multiEncodedC.size() = " << multiEncodedC.size()
           << ", numGroups = " << numGroups << ", numBlocksPerGroup = " << numBlocksPerGroup << endl;

      vector<helib::Ctxt> eDataC;
      eDataC.reserve(numGroups*16);
      for (long j=0; j<numGroups*16; j++) {
        eDataC.push_back(helib::Ctxt(publicKey));
        publicKey.Encrypt(eDataC[j], multiEncodedC[j]);
      }

      // combinedLinearLayer已支持多block, 现在可以直接把numGroups*16个密文
      // 一次性塞给bitLevelRoundTestA, 内部invert/unpack/pack/combinedLinearLayer
      // 都会按eData.size()自动用满线程, 不会有4组并发导致的资源争抢
      double tmC = -NTL::GetTime();
      hAES.bitLevelRoundTestA(eDataC, encryptedAESkeyC);
      tmC += NTL::GetTime();
      cout << "[方案C] 10轮总耗时(" << totalBlocksC << "个block, " << numGroups*16
           << "个密文, 不含最初加密) = " << tmC << " 秒" << endl;

      {
        int correctCount = 0, incorrectCount = 0;
        for (long i=0; i<numGroups*16; i++) {
          if (eDataC[i].isCorrect()) correctCount++;
          else incorrectCount++;
        }
        cout << "[FHE-CORRECTNESS][方案C] isCorrect() 检验: "
             << correctCount << "/" << numGroups*16 << " 通过, "
             << incorrectCount << "/" << numGroups*16 << " 未通过" << endl;
      }

      vector<ZZX> polyC(numGroups*16);
      for (long j=0; j<numGroups*16; j++) secretKey.Decrypt(polyC[j], eDataC[j]);
      Vec<uint8_t> multiCtxtBlocksC(INIT_SIZE, totalBlocksC*16);
      decode4AES_multiGroupBlock16Ctxt(multiCtxtBlocksC.data(), numGroups, numBlocksPerGroup, polyC, ea2);

      unsigned char keyScheduleC[240];
      AESKeyExpansion(keyScheduleC, aesKey.data(), 128);

      long blockMatchCountC = 0, blockMismatchCountC = 0;
      for (long i=0; i<totalBlocksC; i++) {
        uint8_t refBlockC[16];
        Cipher(refBlockC, &multiBlocksC[i*16], keyScheduleC, 10);
        bool blockMatch = true;
        for (int b=0; b<16; b++) {
          if (multiCtxtBlocksC[i*16+b] != refBlockC[b]) { blockMatch = false; break; }
        }
        if (blockMatch) blockMatchCountC++;
        else {
          blockMismatchCountC++;
          if (blockMismatchCountC <= 5)
            cerr << "[方案C] block " << i << " 不匹配!" << endl;
        }
      }
      cout << "[方案C] 逐block正确性: " << blockMatchCountC << "/" << totalBlocksC
           << " 匹配, " << blockMismatchCountC << "/" << totalBlocksC << " 不匹配" << endl;
      if (blockMismatchCountC == 0)
        cout << "[方案C] 结果: 全部一致 ✓" << endl;
      else
        cerr << "[方案C] 结果: 存在不一致 ✗" << endl;
    }
    // ===== [方案C] 测试结束 =====
    cout << "\n\n######## [MARKER] >>>>> 方案C 结束 <<<<< ########\n" << endl;
    cerr << "\n\n######## [MARKER] >>>>> 方案C 结束 <<<<< ########\n" << endl;

    // ===== [NOISE-STATS] 方案3：蒙特卡洛统计解密错误率(密钥+明文都随机) =====
    // 注意：这里跟boot=0/boot=1共用同一套流程，boot=1时每次trial都会重新走一遍
    // 完整的自举同态AES流程，单次耗时会比boot=0慢很多，trialCount先设小一点验证。
    {
      cerr << "\n===== [NOISE-STATS] 方案3: 蒙特卡洛统计解密错误率(密钥+明文都随机) =====" << endl;
      long trialCount = 1; // 临时改小以快速验证方案B，之后再改回50正式统计
      long totalBitErrors = 0;
      long totalBits = 0;
      long totalTrialsWithError = 0;
      double totalTimeSec = 0.0;
      long trialNumBlocks = ea2.size();
      for (long trial = 0; trial < trialCount; trial++) {
        cerr << "[NOISE-STATS] --- 开始第 " << (trial+1) << "/" << trialCount
             << " 次独立试验 (" << trialNumBlocks << "-block 批量) ---" << endl;
        double tm = -NTL::GetTime();
        Vec<uint8_t> trialAesKey(INIT_SIZE, 16);
        Vec<uint8_t> trialBlocks(INIT_SIZE, trialNumBlocks*16);
        {
          GF2X rnd;
          random(rnd, 8 * trialAesKey.length());
          BytesFromGF2X(trialAesKey.data(), rnd, trialAesKey.length());
          GF2X rnd2;
          random(rnd2, 8 * trialBlocks.length());
          BytesFromGF2X(trialBlocks.data(), rnd2, trialBlocks.length());
        }
        vector<vector<helib::Ctxt>> trialEncryptedAESkeyB;
        hAES.encryptAESkey_broadcastBlock16Ctxt(trialEncryptedAESkeyB, trialAesKey, publicKey);
        vector<ZZX> trialEncoded;
        encode4AES_multiBlock16Ctxt(trialEncoded, trialBlocks.data(), trialNumBlocks, ea2);
        vector<helib::Ctxt> trialEData;
        trialEData.reserve(16);
        for (long j = 0; j < 16; j++) {
          trialEData.push_back(helib::Ctxt(publicKey));
          publicKey.Encrypt(trialEData[j], trialEncoded[j]);
        }
        hAES.bitLevelRoundTestA(trialEData, trialEncryptedAESkeyB);
        tm += NTL::GetTime();
        totalTimeSec += tm;
        cerr << "[NOISE-STATS] 第" << (trial+1) << "次试验耗时(含密钥重新加密, "
             << trialNumBlocks << "个block) = " << tm << " 秒" << endl;
        vector<ZZX> trialPoly(16);
        for (long j = 0; j < 16; j++) secretKey.Decrypt(trialPoly[j], trialEData[j]);
        Vec<uint8_t> trialCtxtBlocks(INIT_SIZE, trialNumBlocks*16);
        decode4AES_multiBlock16Ctxt(trialCtxtBlocks.data(), trialNumBlocks, trialPoly, hAES.getEA());
        unsigned char trialKeySchedule[240];
        AESKeyExpansion(trialKeySchedule, trialAesKey.data(), 128);
        bool trialHasError = false;
        long trialBlockMismatch = 0;
        for (long i = 0; i < trialNumBlocks; i++) {
          uint8_t trialRefBlock[16];
          Cipher(trialRefBlock, &trialBlocks[i*16], trialKeySchedule, 10);
          bool blockHasError = false;
          for (long b = 0; b < 16; b++) {
            uint8_t xorVal = trialCtxtBlocks[i*16+b] ^ trialRefBlock[b];
            for (long bit = 0; bit < 8; bit++) {
              totalBits++;
              if ((xorVal >> bit) & 1) {
                totalBitErrors++;
                blockHasError = true;
              }
            }
          }
          if (blockHasError) { trialHasError = true; trialBlockMismatch++; }
        }
        if (trialHasError) {
          totalTrialsWithError++;
          cerr << "[NOISE-STATS] 第" << (trial+1) << "次试验: 出现解密错误! ✗ ("
               << trialBlockMismatch << "/" << trialNumBlocks << " block不匹配)" << endl;
        } else {
          cerr << "[NOISE-STATS] 第" << (trial+1) << "次试验: 结果正确 ✓ ("
               << trialNumBlocks << "/" << trialNumBlocks << " block全部匹配)" << endl;
        }
      }


      double bitErrorRate = (double)totalBitErrors / (double)totalBits;
      double trialErrorRate = (double)totalTrialsWithError / (double)trialCount;
      cerr << "[NOISE-STATS] ===== 方案3 汇总 =====" << endl;
      cerr << "[NOISE-STATS] 蒙特卡洛试验次数 = " << trialCount << endl;
      cerr << "[NOISE-STATS] 总耗时 = " << totalTimeSec << " 秒" << endl;
      cerr << "[NOISE-STATS] 总bit数 = " << totalBits
           << ", 出错bit数 = " << totalBitErrors << endl;
      cerr << "[NOISE-STATS] 逐bit错误率 = " << bitErrorRate << endl;
      cerr << "[NOISE-STATS] 逐trial(整块)错误率 = " << trialErrorRate << endl;
    }

  }

  return 0;
}
