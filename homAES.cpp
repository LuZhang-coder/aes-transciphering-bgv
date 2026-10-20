/** homAES.cpp - adapted for HElib 2.2.0 */
#include <cstdlib>
#include <thread>
#include <atomic>
#include <chrono>
#include <cstring>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/mat_GF2E.h>
#include <NTL/BasicThreadPool.h>
#include "homAES.h"
#include "combined_matrix.h"
#include "decrypt_matrix_v2.h"

using namespace std;
extern int getSBoxValue(int num);
uint8_t gf256mul(uint8_t a, uint8_t b)
{
  uint8_t p = 0;
  for (int i = 0; i < 8; i++) {
    if (b & 1) p ^= a;
    uint8_t hiBitSet = a & 0x80;
    a <<= 1;
    if (hiBitSet) a ^= 0x1B;
    b >>= 1;
  }
  return p;
}
// GF(2^8)乘法逆元(纯求逆, 不含affine), AES标准域(不可约多项式0x11B)
// getInverseValue(0)=0 (AES约定)
int getInverseValue(int x)
{
  uint8_t xb = (uint8_t)x;
  if (xb == 0) return 0;
  uint8_t result = 1;
  uint8_t base = xb;
  int exp = 254;
  while (exp > 0) {
    if (exp & 1) result = gf256mul(result, base);
    base = gf256mul(base, base);
    exp >>= 1;
  }
  return (int)result;
}
using namespace NTL;

static void buildAffine(vector<PolyType>& binMat, PolyType* binVec,
                        const unsigned char cc[],
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);

inline void buildAffineEnc(vector<PolyType>& binMat, PolyType& binVec,
                           const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  static unsigned char cc[] = { 31, 62, 124, 248, 241, 227, 199, 143, 99 };
  buildAffine(binMat, &binVec, cc, ea2);
}

inline void buildAffineDec(vector<PolyType>& binMat,
                           const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  static unsigned char cc[] = { 74, 148, 41, 82, 164, 73, 146, 37, 0 };
  buildAffine(binMat, nullptr, cc, ea2);
}

static void buildLinEnc(vector<PolyType>& encLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void buildLinDec(vector<PolyType>& decLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void encRowColTran(helib::Ctxt& c, const vector<PolyType>& encLinTran,
                          const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void encRowShift(helib::Ctxt& c, const vector<PolyType>& encLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void decRowColTran(helib::Ctxt& c, const vector<PolyType>& decLinTran,
                          const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void decRowShift(helib::Ctxt& c, const vector<PolyType>& decLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void invert(vector<helib::Ctxt>& data);

// ===== 方案A新增：byte<->bit 验证用辅助函数 =====
static void buildByteToBitMat(vector<PolyType>& bitMat,
                               const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void unpackByteToBits(vector<helib::Ctxt>& bits, const helib::Ctxt& byteCtxt,
                              const vector<PolyType>& bitMat,
                              const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void packBitsToByte(helib::Ctxt& byteCtxt, const vector<helib::Ctxt>& bits,
                            const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
static void bitLevelShiftRows(vector<vector<helib::Ctxt>>& stateBits);
static void bitLevelMixColumns(vector<vector<helib::Ctxt>>& stateBits);
// ===== 方案A新增结束 =====
static void packCtxt(vector<helib::Ctxt>& to, const vector<helib::Ctxt>& from,
                     const GF2X& XinSlots);
static void unackCtxt(vector<helib::Ctxt>& to, const vector<helib::Ctxt>& from,
                      const Mat<GF2X>& unpackConsts);

static const uint8_t aesPolyBytes[] = { 0x1B, 0x1 };
const GF2X HomAES::aesPoly = GF2XFromBytes(aesPolyBytes, 2);

HomAES::HomAES(const helib::Context& context)
  : ea2(context, aesPoly, context.getAlMod())
#ifndef USE_ZZX_POLY
  , affVec(0, context, helib::IndexSet(0, context.numPrimes()-1))
#endif
{
  helib::assertEq(context.getZMStar().OrderOf(0) % 16, 0l,
                  "First dimension order must be divisible by 16");
  buildAffineEnc(encAffMat, affVec, ea2);
  buildAffineDec(decAffMat, ea2);
  buildLinEnc(encLinTran, ea2);
  buildLinDec(decLinTran, ea2);
  if (context.isBootstrappable())
    setPackingConstants();
}

void HomAES::encryptAESkey(vector<helib::Ctxt>& eKey, Vec<uint8_t>& aesKey,
                           const helib::PubKey& hePK) const
{
  helib::assertTrue<helib::InvalidArgument>(
    aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
    "Key length must be 16, 24, or 32");
  extern long AESKeyExpansion(unsigned char[], unsigned char[], int);
  uint8_t roundKeySchedule[240];
  long nRoundKeys = AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);
  long blocksPerCtxt = ea2.size() / 16;
  Vec<uint8_t> expanded(INIT_SIZE, nRoundKeys*blocksPerCtxt*16);
  for (long i=0; i<nRoundKeys; i++) {
    uint8_t* roundKey = &roundKeySchedule[16*i];
    for (long j=0; j<blocksPerCtxt; j++)
      memcpy(&expanded[16*(i*blocksPerCtxt+j)], roundKey, 16);
  }
  Vec<ZZX> encoded;
  encode4AES(encoded, expanded, ea2);
  { helib::Ctxt tmpCtxt(hePK);
    eKey.resize(encoded.length(), tmpCtxt); }
  for (long i=0; i<(long)eKey.size(); i++)
    hePK.Encrypt(eKey[i], encoded[i]);
}
void HomAES::encryptAESkey_singleBlock16Ctxt(vector<vector<helib::Ctxt>>& eKeySchedule,
                                                Vec<uint8_t>& aesKey,
                                                const helib::PubKey& hePK) const
{
  helib::assertTrue<helib::InvalidArgument>(
    aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
    "Key length must be 16, 24, or 32");
  extern long AESKeyExpansion(unsigned char[], unsigned char[], int);
  uint8_t roundKeySchedule[240];
  long nRoundKeys = AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);
  eKeySchedule.resize(nRoundKeys);
  for (long r=0; r<nRoundKeys; r++) {
    vector<ZZX> encoded16;
    encode4AES_singleBlock16Ctxt(encoded16, &roundKeySchedule[16*r], ea2);
    helib::Ctxt tmpCtxt(hePK);
    eKeySchedule[r].resize(16, tmpCtxt);
    for (long k=0; k<16; k++)
      hePK.Encrypt(eKeySchedule[r][k], encoded16[k]);
  }
}
void HomAES::encryptAESkey_broadcastBlock16Ctxt(vector<vector<helib::Ctxt>>& eKeySchedule,
                                                Vec<uint8_t>& aesKey,
                                                const helib::PubKey& hePK) const
{
  helib::assertTrue<helib::InvalidArgument>(
    aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
    "Key length must be 16, 24, or 32");
  extern long AESKeyExpansion(unsigned char[], unsigned char[], int);
  uint8_t roundKeySchedule[240];
  long nRoundKeys = AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);
  eKeySchedule.resize(nRoundKeys);
  for (long r=0; r<nRoundKeys; r++) {
    vector<ZZX> encoded16;
    encode4AES_broadcastBlock16Ctxt(encoded16, &roundKeySchedule[16*r], ea2);
    helib::Ctxt tmpCtxt(hePK);
    eKeySchedule[r].resize(16, tmpCtxt);
    for (long k=0; k<16; k++)
      hePK.Encrypt(eKeySchedule[r][k], encoded16[k]);
  }
}
void HomAES::encryptAESkey_multiGroupBroadcastCtxt(vector<vector<helib::Ctxt>>& eKeySchedule,
                                                Vec<uint8_t>& aesKey, long numGroups,
                                                const helib::PubKey& hePK) const
{
  helib::assertTrue<helib::InvalidArgument>(
    aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
    "Key length must be 16, 24, or 32");
  extern long AESKeyExpansion(unsigned char[], unsigned char[], int);
  uint8_t roundKeySchedule[240];
  long nRoundKeys = AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);
  eKeySchedule.resize(nRoundKeys);
  for (long r=0; r<nRoundKeys; r++) {
    vector<ZZX> encoded16;
    encode4AES_broadcastBlock16Ctxt(encoded16, &roundKeySchedule[16*r], ea2);
    helib::Ctxt tmpCtxt(hePK);
    eKeySchedule[r].resize(numGroups*16, tmpCtxt);
    for (long g=0; g<numGroups; g++)
      for (long k=0; k<16; k++)
        hePK.Encrypt(eKeySchedule[r][g*16 + k], encoded16[k]);
  }
}


void HomAES::encryptAESkey_bitLevel(vector<vector<vector<helib::Ctxt>>>& eKeyBits,
                                     Vec<uint8_t>& aesKey,
                                     const helib::PubKey& hePK) const
{
  helib::assertTrue<helib::InvalidArgument>(
    aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
    "Key length must be 16, 24, or 32");
  extern long AESKeyExpansion(unsigned char[], unsigned char[], int);
  uint8_t roundKeySchedule[240];
  long nRoundKeys = AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);
  long numSlots = ea2.size();
  eKeyBits.resize(nRoundKeys);
  for (long r=0; r<nRoundKeys; r++) {
    eKeyBits[r].resize(16);
    for (long j=0; j<16; j++) {
      eKeyBits[r][j].reserve(8);
      uint8_t byteVal = roundKeySchedule[16*r + j];
      for (int k=0; k<8; k++) {
        uint8_t bitVal = (byteVal >> k) & 1;
        Vec<uint8_t> bitVec(INIT_SIZE, numSlots);
        for (long s=0; s<numSlots; s++) bitVec[s] = bitVal;
        Vec<ZZX> bitEncoded;
        encode4AES(bitEncoded, bitVec, ea2);
        helib::Ctxt bitCtxt(hePK);
        hePK.Encrypt(bitCtxt, bitEncoded[0]);
        eKeyBits[r][j].push_back(bitCtxt);
      }
    }
  }
}
void HomAES::setPackingConstants()
{
  const helib::Context& context = ea2.getContext();
  const helib::EncryptedArrayDerived<helib::PA_GF2>& ea =
    context.getEA().getDerived(helib::PA_GF2());
  { vector<GF2X> slots(ea.size(), GF2X(1,1));
    ZZX tmp; ea.encode(tmp, slots);
    conv(XinSlots, tmp); }
  long e = ea.getDegree() / 8;
  std::cerr << "[诊断] 自然ea.getDegree()=" << ea.getDegree() << " e(=ea.getDegree()/8)=" << e << std::endl;
  GF2EBak bak; bak.save();
  GF2XModulus F0(ea.getTab().getFactors()[0]);
  GF2E::init(F0);
  mat_GF2E Kinv(INIT_SIZE, e, e);
  for (long j=0; j<e; j++) conv(Kinv[0][j], GF2X(j,1));
  for (long i=1; i<e; i++) for (long j=0; j<e; j++)
    power(Kinv[i][j], Kinv[0][j], 1L<<(8*i));
  mat_GF2E K(INIT_SIZE, e, e);
  inv(K, Kinv);
  unpacking.SetDims(e,e);
  for (long i=0; i<e; i++) for (long j=0; j<e; j++) {
    vector<GF2X> slots(ea.size(), rep(K[i][j]));
    ZZX tmp; ea.encode(tmp, slots);
    conv(unpacking[i][j], tmp);
  }
}

void HomAES::homAESenc(vector<helib::Ctxt>& eData, const vector<helib::Ctxt>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return;
  for (long j=0; j<(long)eData.size(); j++) eData[j] += aesKey[0];
  double t_invert=0, t_affine=0, t_shift=0, t_boot=0;
  double wall_boot=0, wall_invert=0, wall_unpack=0, wall_linear=0, wall_pack=0;
  std::chrono::steady_clock::time_point _wnow0;
  for (long i=1; i<(long)aesKey.size(); i++) {
    double t0, t1;
    std::cerr << "[ENC pre-invert] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 80) { auto _bw0=std::chrono::steady_clock::now(); t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); batchRecrypt(eData); t_boot+=NTL::GetTime()-t0; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count(); }
    auto _wall0 = std::chrono::steady_clock::now();
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); invert(eData); t1=NTL::GetTime(); t_invert+=t1-t0;
    auto _wall1 = std::chrono::steady_clock::now();
    double _wallSec = std::chrono::duration<double>(_wall1-_wall0).count();
    wall_invert += _wallSec;
    std::cerr << "  [WALL] invert_wall=" << _wallSec << "s (NTL_time=" << (t1-t0) << "s)" << std::endl;
    std::cerr << "  [TIME] round=" << i << " invert=" << t1-t0 << "s" << std::endl;
    std::cerr << "[ENC pre-affine] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 35) { auto _bw0=std::chrono::steady_clock::now(); t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); batchRecrypt(eData); t_boot+=NTL::GetTime()-t0; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count(); }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    for (long j=0; j<(long)eData.size(); j++) {
      helib::applyLinPolyLL(eData[j], encAffMat, ea2.getDegree());
      eData[j].addConstant(affVec);
    }

    t1=NTL::GetTime(); t_affine+=t1-t0;
    if (i == 1) {
      std::cerr << "  [pack/unpack test] before unpack bitCap=" << eData[0].bitCapacity() << std::endl;
      double t_unpack_test = NTL::GetTime();
      vector<helib::Ctxt> unpacked;
      unackCtxt(unpacked, eData, unpacking);
      double unpack_elapsed = NTL::GetTime() - t_unpack_test;
      std::cerr << "  [unpack test] elapsed=" << unpack_elapsed
                 << "s, bitCap[0]=" << unpacked[0].bitCapacity()
                 << " count=" << unpacked.size() << std::endl;
      double t_pack_test = NTL::GetTime();
      vector<helib::Ctxt> repacked;
      packCtxt(repacked, unpacked, XinSlots);
      double pack_elapsed = NTL::GetTime() - t_pack_test;
      std::cerr << "  [pack test] elapsed=" << pack_elapsed
                 << "s, bitCap[0]=" << repacked[0].bitCapacity()
                 << " count=" << repacked.size() << std::endl;
    }

    std::cerr << "  [TIME] round=" << i << " affine=" << t1-t0 << "s" << std::endl;
    std::cerr << "[ENC pre-shift] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 35) { auto _bw0=std::chrono::steady_clock::now(); t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); batchRecrypt(eData); t_boot+=NTL::GetTime()-t0; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count(); }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    if (i<(long)aesKey.size()-1)
      for (long j=0; j<(long)eData.size(); j++) encRowColTran(eData[j], encLinTran, ea2);
    else
      for (long j=0; j<(long)eData.size(); j++) encRowShift(eData[j], encLinTran, ea2);
    for (long j=0; j<(long)eData.size(); j++) eData[j] += aesKey[i];
    t1=NTL::GetTime(); t_shift+=t1-t0;
    std::cerr << "  [TIME] round=" << i << " shift+keyadd=" << t1-t0 << "s" << std::endl;
  }
  std::cerr << "[ENC TOTAL] invert=" << t_invert << "s  affine=" << t_affine << "s  shift=" << t_shift << "s  boot=" << t_boot << "s" << std::endl;
}

void HomAES::homAESenc(vector<helib::Ctxt>& eData, const vector<helib::Ctxt>& aesKey,
                       const Vec<uint8_t> inBytes) const
{
  Vec<ZZX> encodedBytes;
  encode4AES(encodedBytes, inBytes, ea2);
  eData.resize(encodedBytes.length(), helib::Ctxt(helib::ZeroCtxtLike, aesKey[0]));
  for (long i=0; i<(long)eData.size(); i++) eData[i].DummyEncrypt(encodedBytes[i]);
  homAESenc(eData, aesKey);
}

void HomAES::homAESdec(vector<helib::Ctxt>& eData, const vector<helib::Ctxt>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return;
  for (long i=aesKey.size()-1; i>0; i--) {
    for (long j=0; j<(long)eData.size(); j++) eData[j] -= aesKey[i];
    std::cerr << "[ENC pre-shift] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 35) batchRecrypt(eData);
    if (i<(long)aesKey.size()-1)
      for (long j=0; j<(long)eData.size(); j++) decRowColTran(eData[j], decLinTran, ea2);
    else
      for (long j=0; j<(long)eData.size(); j++) decRowShift(eData[j], decLinTran, ea2);
    if (eData[0].bitCapacity() < 35) batchRecrypt(eData);
    for (long j=0; j<(long)eData.size(); j++) {
      eData[j].addConstant(affVec);
      helib::applyLinPolyLL(eData[j], decAffMat, ea2.getDegree());
    }
    std::cerr << "[ENC pre-invert] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 80) batchRecrypt(eData);
    invert(eData);
  }
  for (long j=0; j<(long)eData.size(); j++) eData[j] -= aesKey[0];
}

void HomAES::homAESdec(vector<helib::Ctxt>& eData, const vector<helib::Ctxt>& aesKey,
                       const Vec<uint8_t> inBytes) const
{
  Vec<ZZX> encodedBytes;
  encode4AES(encodedBytes, inBytes, ea2);
  eData.resize(encodedBytes.length(), helib::Ctxt(helib::ZeroCtxtLike, aesKey[0]));
  for (long i=0; i<(long)eData.size(); i++) eData[i].DummyEncrypt(encodedBytes[i]);
  homAESdec(eData, aesKey);
}

void HomAES::batchRecrypt(vector<helib::Ctxt>& data) const
{
  HELIB_TIMER_START;
  const helib::PubKey& pk = data[0].getPubKey();
  if (!pk.isBootstrappable()) return;
  if (data.size()>1 && unpacking.NumRows()==0) return;
  vector<helib::Ctxt>* pData = &data;
  vector<helib::Ctxt> fullyPacked;
  if (data.size()>1) {
    std::cerr << "  [pack test] before, bitCap[0]=" << data[0].bitCapacity() << " count=" << data.size() << std::endl;
    double t_pack = NTL::GetTime();
    packCtxt(fullyPacked, data, XinSlots);
    std::cerr << "  [pack test] elapsed=" << NTL::GetTime()-t_pack
               << "s, after bitCap[0]=" << fullyPacked[0].bitCapacity()
               << " count=" << fullyPacked.size() << std::endl;
    pData = &fullyPacked;
  }
  HELIB_NTIMER_START(recryption);
  auto _rBody = [&](size_t j) {
#ifdef USE_PAPER_VERSION
    pk.reCrypt((*pData)[j], true);
#else
    pk.reCrypt((*pData)[j]);
#endif
  };
  if (pData->size() < 8) {
    for (size_t j=0; j<pData->size(); j++) _rBody(j);
  } else {
    // 分批并发: 每批用 _rBatchSize 个std::thread并发调用reCrypt, 批与批之间顺序执行,
    // 避免外层并发数过多时跟NTL内部线程池产生过度资源竞争。
    // 可通过环境变量 RECRYPT_BATCH_SIZE 调整批大小(不设置则默认全量并发, 等价于原来的行为)。
    size_t _rBatchSize = pData->size();
    const char* _envBatch = std::getenv("RECRYPT_BATCH_SIZE");
    if (_envBatch != nullptr) {
      long _v = std::atol(_envBatch);
      if (_v > 0) _rBatchSize = (size_t)_v;
    }
    std::cerr << "[RECRYPT-TUNE] pData->size()=" << pData->size()
              << " batchSize=" << _rBatchSize << std::endl;
    for (size_t base = 0; base < pData->size(); base += _rBatchSize) {
      size_t end = std::min(base + _rBatchSize, pData->size());
      std::vector<std::thread> _rThreads;
      for (size_t j = base; j < end; j++) _rThreads.emplace_back(_rBody, j);
      for (auto& t : _rThreads) t.join();
    }
  }
  HELIB_NTIMER_STOP(recryption);
  if (fullyPacked.size()>0) {
    std::cerr << "  [unpack test] before, bitCap[0]=" << fullyPacked[0].bitCapacity() << " count=" << fullyPacked.size() << std::endl;
    double t_unpack = NTL::GetTime();
    unackCtxt(data, fullyPacked, unpacking);
    std::cerr << "  [unpack test] elapsed=" << NTL::GetTime()-t_unpack
               << "s, after bitCap[0]=" << data[0].bitCapacity()
               << " count=" << data.size() << std::endl;
  }
}

static void buildAffine(vector<PolyType>& binMat, PolyType* binVec,
                        const unsigned char cc[],
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  vector<GF2X> scratch(8);
  for (long j=0; j<8; j++) GF2XFromBytes(scratch[j], &cc[j], 1);
  vector<GF2X> C;
  ea2.buildLinPolyCoeffs(C, scratch);
  vector<ZZX> zzxMat(8);
  scratch.resize(ea2.size());
  for (long j=0; j<8; j++) {
    for (long i=0; i<ea2.size(); i++) scratch[i] = C[j];
    ea2.encode(zzxMat[j], scratch);
  }
#ifndef USE_ZZX_POLY
  binMat.resize(8, helib::DoubleCRT(0, ea2.getContext(),
                helib::IndexSet(0, ea2.getContext().numPrimes()-1)));
  for (long j=0; j<8; j++) binMat[j] = zzxMat[j];
#else
  binMat = zzxMat;
#endif
  if (binVec != nullptr) {
    GF2X cc8; GF2XFromBytes(cc8, &cc[8], 1);
    for (long i=0; i<ea2.size(); i++) scratch[i] = cc8;
    ZZX tmpZZX; ea2.encode(tmpZZX, scratch);
    *binVec = tmpZZX;
  }
}

static void buildLinEnc(vector<PolyType>& encLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  uint8_t hexVals[3] = { 1, 2, 3 };
  Vec<uint8_t> bytes(INIT_SIZE, ea2.size());
  memset(bytes.data(), 0, bytes.length());
  long blocksPerCtxt = ea2.size() / 16;
  Vec<ZZX> tmp;
  encLinTran.resize(6, helib::DoubleCRT(0, ea2.getContext(),
                    helib::IndexSet(0, ea2.getContext().numPrimes()-1)));
  for (long i=0; i<3; i++) {
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j];
      bptr[0]=bptr[4]=bptr[8]=bptr[12]=hexVals[i];
    }
    encode4AES(tmp, bytes, ea2); encLinTran[i] = tmp[0];
  }
  for (long i=0; i<3; i++) {
    memset(bytes.data(), 0, bytes.length());
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j+i+1];
      bptr[0]=bptr[4]=bptr[8]=bptr[12]=1;
    }
    encode4AES(tmp, bytes, ea2); encLinTran[i+3] = tmp[0];
  }
}

static void encRowColTran(helib::Ctxt& c, const vector<PolyType>& encLinTran,
                          const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long rotAmount = ea2.getContext().getZMStar().OrderOf(0) / 16;
  c.cleanUp();
  std::cerr << "  [mixcol] entry bitCap=" << c.bitCapacity() << std::endl;
  helib::Ctxt c1(c), c6(c), c11(c);
  ea2.rotate1D(c1,0,rotAmount); ea2.rotate1D(c6,0,6*rotAmount); ea2.rotate1D(c11,0,11*rotAmount);
  c1.cleanUp(); c6.cleanUp(); c11.cleanUp();
  std::cerr << "  [mixcol] after rotate bitCap=" << c1.bitCapacity() << std::endl;
  const PolyType& p1=encLinTran[0]; const PolyType& p2=encLinTran[1]; const PolyType& p3=encLinTran[2];
  helib::Ctxt sum(c), tmpA(c1), tmpB(c11);
  tmpA += c6; sum.multByConstant(p2); tmpA.multByConstant(p1); tmpB.multByConstant(p3);
  std::cerr << "  [mixcol] after first multByConstant bitCap=" << sum.bitCapacity() << std::endl;
  sum += tmpA; sum += tmpB;
  tmpA=c; tmpA+=c1; tmpA.multByConstant(p1);
  tmpB=c6; tmpB.multByConstant(p3);
  helib::Ctxt tmpC(c11); tmpC.multByConstant(p2);
  tmpA+=tmpB; tmpA+=tmpC; ea2.rotate1D(tmpA,0,rotAmount); sum+=tmpA;
  tmpA=c; tmpA+=c11; tmpA.multByConstant(p1);
  tmpB=c1; tmpB.multByConstant(p3); tmpC=c6; tmpC.multByConstant(p2);
  tmpA+=tmpB; tmpA+=tmpC; ea2.rotate1D(tmpA,0,2*rotAmount); sum+=tmpA;
  c.multByConstant(p3); c1.multByConstant(p2); c6+=c11; c6.multByConstant(p1);
  c+=c1; c+=c6; ea2.rotate1D(c,0,3*rotAmount); c+=sum; c.cleanUp();
  std::cerr << "  [mixcol] exit bitCap=" << c.bitCapacity() << std::endl;
}
static void encRowShift(helib::Ctxt& c, const vector<PolyType>& encLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long rotAmount = ea2.getContext().getZMStar().OrderOf(0) / 16;
  c.cleanUp();
  helib::Ctxt c4(c), c8(c), c12(c);
  ea2.rotate1D(c4,0,4*rotAmount); ea2.rotate1D(c8,0,8*rotAmount); ea2.rotate1D(c12,0,12*rotAmount);
  c4.cleanUp(); c8.cleanUp(); c12.cleanUp();
  const PolyType& p0=encLinTran[0]; const PolyType& p4=encLinTran[3];
  const PolyType& p8=encLinTran[4]; const PolyType& p12=encLinTran[5];
  c.multByConstant(p0); c4.multByConstant(p12); c+=c4;
  c8.multByConstant(p8); c+=c8; c12.multByConstant(p4); c+=c12;
}

static void buildLinDec(vector<PolyType>& decLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  uint8_t hexVals[4] = { 0x9, 0xB, 0xD, 0xE };
  Vec<uint8_t> bytes(INIT_SIZE, ea2.size());
  memset(bytes.data(), 0, bytes.length());
  long blocksPerCtxt = ea2.size() / 16;
  Vec<ZZX> tmp;
  decLinTran.resize(8, helib::DoubleCRT(0, ea2.getContext(),
                    helib::IndexSet(0, ea2.getContext().numPrimes()-1)));
  for (long i=0; i<4; i++) {
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j];
      bptr[0]=bptr[4]=bptr[8]=bptr[12]=hexVals[i];
    }
    encode4AES(tmp, bytes, ea2); decLinTran[i] = tmp[0];
  }
  for (long i=0; i<4; i++) {
    memset(bytes.data(), 0, bytes.length());
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j+i];
      bptr[0]=bptr[4]=bptr[8]=bptr[12]=1;
    }
    encode4AES(tmp, bytes, ea2); decLinTran[i+4] = tmp[0];
  }
}

static void decRowColTran(helib::Ctxt& c, const vector<PolyType>& decLinTran,
                          const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long rotAmount = ea2.getContext().getZMStar().OrderOf(0) / 16;
  c.cleanUp();
  helib::Ctxt cf(c), ce(c), cd(c);
  ea2.rotate1D(cf,0,15*rotAmount); ea2.rotate1D(ce,0,14*rotAmount); ea2.rotate1D(cd,0,13*rotAmount);
  cf.cleanUp(); ce.cleanUp(); cd.cleanUp();
  const PolyType& p9=decLinTran[0]; const PolyType& pB=decLinTran[1];
  const PolyType& pD=decLinTran[2]; const PolyType& pE=decLinTran[3];
  helib::Ctxt sum(c), tmpF(cf), tmpE(ce), tmpD(cd);
  sum.multByConstant(pE); tmpF.multByConstant(pB); tmpE.multByConstant(pD); tmpD.multByConstant(p9);
  sum+=tmpF; sum+=tmpE; sum+=tmpD;
  helib::Ctxt tmp(c); tmpF=cf; tmpE=ce; tmpD=cd;
  tmp.multByConstant(p9); tmpF.multByConstant(pE); tmpE.multByConstant(pB); tmpD.multByConstant(pD);
  tmp+=tmpF; tmp+=tmpE; tmp+=tmpD; ea2.rotate1D(tmp,0,5*rotAmount); sum+=tmp;
  tmp=c; tmpF=cf; tmpE=ce; tmpD=cd;
  tmp.multByConstant(pD); tmpF.multByConstant(p9); tmpE.multByConstant(pE); tmpD.multByConstant(pB);
  tmp+=tmpF; tmp+=tmpE; tmp+=tmpD; ea2.rotate1D(tmp,0,10*rotAmount); sum+=tmp;
  c.multByConstant(pB); cf.multByConstant(pD); ce.multByConstant(p9); cd.multByConstant(pE);
  c+=cf; c+=ce; c+=cd; ea2.rotate1D(c,0,15*rotAmount); c+=sum; c.cleanUp();
}

static void decRowShift(helib::Ctxt& c, const vector<PolyType>& decLinTran,
                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long rotAmount = ea2.getContext().getZMStar().OrderOf(0) / 16;
  c.cleanUp();
  helib::Ctxt c4(c), c8(c), c12(c);
  ea2.rotate1D(c4,0,4*rotAmount); ea2.rotate1D(c8,0,8*rotAmount); ea2.rotate1D(c12,0,12*rotAmount);
  c4.cleanUp(); c8.cleanUp(); c12.cleanUp();
  const PolyType& p0=decLinTran[4]; const PolyType& p4=decLinTran[5];
  const PolyType& p8=decLinTran[6]; const PolyType& p12=decLinTran[7];
  c.multByConstant(p0); c4.multByConstant(p4); c+=c4;
  c8.multByConstant(p8); c+=c8; c12.multByConstant(p12); c+=c12;
}

// ===== 方案A 单block验证专用编码/解码 =====
// 把1个AES block(16字节)编码成16个密文(每个密文对应一个字节位置)，
// 每个密文只有slot 0放真实字节值，其余959个slot填0(占位，保证密文结构完整)
void encode4AES_singleBlock16Ctxt(vector<ZZX>& encData, const uint8_t block[16],
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  encData.resize(16);
  for (long k=0; k<16; k++) {
    vector<GF2X> slots(ea2.size(), GF2X::zero());
    GF2XFromBytes(slots[0], &block[k], 1);
    ea2.encode(encData[k], slots);
  }
}
// 对应解码: 从16个密文里各取slot 0，还原成1个block的16字节
void decode4AES_singleBlock16Ctxt(uint8_t block[16], const vector<ZZX>& encData,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  vector<GF2X> slots;
  for (long k=0; k<16; k++) {
    ea2.decode(slots, encData[k]);
    BytesFromGF2X(&block[k], slots[0], 1);
  }
}
void encode4AES_broadcastBlock16Ctxt(vector<ZZX>& encData, const uint8_t block[16],
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long numSlots = ea2.size();
  encData.resize(16);
  for (long k=0; k<16; k++) {
    GF2X byteVal;
    GF2XFromBytes(byteVal, &block[k], 1);
    vector<GF2X> slots(numSlots, byteVal);
    ea2.encode(encData[k], slots);
  }
}
void encode4AES_multiGroupBlock16Ctxt(vector<ZZX>& encData, const uint8_t* blocks,
                long numGroups, long numBlocksPerGroup,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long numSlots = ea2.size();
  helib::assertTrue<helib::InvalidArgument>(numBlocksPerGroup <= numSlots,
    "numBlocksPerGroup must not exceed numSlots");
  encData.resize(numGroups * 16);
  for (long g=0; g<numGroups; g++) {
    const uint8_t* groupBlocks = blocks + (size_t)g * numBlocksPerGroup * 16;
    for (long k=0; k<16; k++) {
      vector<GF2X> slots(numSlots, GF2X::zero());
      for (long i=0; i<numBlocksPerGroup; i++) {
        GF2XFromBytes(slots[i], &groupBlocks[i*16 + k], 1);
      }
      ea2.encode(encData[g*16 + k], slots);
    }
  }
}
void decode4AES_multiGroupBlock16Ctxt(uint8_t* blocks, long numGroups, long numBlocksPerGroup,
                const vector<ZZX>& encData,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  vector<GF2X> slots;
  for (long g=0; g<numGroups; g++) {
    uint8_t* groupBlocks = blocks + (size_t)g * numBlocksPerGroup * 16;
    for (long k=0; k<16; k++) {
      ea2.decode(slots, encData[g*16 + k]);
      for (long i=0; i<numBlocksPerGroup; i++) {
        BytesFromGF2X(&groupBlocks[i*16 + k], slots[i], 1);
      }
    }
  }
}
void encode4AES_multiBlock16Ctxt(vector<ZZX>& encData, const uint8_t* blocks, long numBlocks,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long numSlots = ea2.size();
  helib::assertTrue<helib::InvalidArgument>(numBlocks <= numSlots,
    "numBlocks must not exceed numSlots");
  encData.resize(16);
  for (long k=0; k<16; k++) {
    vector<GF2X> slots(numSlots, GF2X::zero());
    for (long i=0; i<numBlocks; i++) {
      GF2XFromBytes(slots[i], &blocks[i*16 + k], 1);
    }
    ea2.encode(encData[k], slots);
  }
}
void decode4AES_multiBlock16Ctxt(uint8_t* blocks, long numBlocks, const vector<ZZX>& encData,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  vector<GF2X> slots;
  for (long k=0; k<16; k++) {
    ea2.decode(slots, encData[k]);
    for (long i=0; i<numBlocks; i++) {
      BytesFromGF2X(&blocks[i*16 + k], slots[i], 1);
    }
  }
}
void encode4AES(Vec<ZZX>& encData, const Vec<uint8_t>& data,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long nBlocks = helib::divc(data.length(), 16);
  long blocksPerCtxt = ea2.size() / 16;
  long nCtxt = helib::divc(nBlocks, blocksPerCtxt);
  encData.SetLength(nCtxt);
  for (long i=0; i<nCtxt; i++) {
    vector<GF2X> slots(ea2.size(), GF2X::zero());
    for (long j=0; j<blocksPerCtxt; j++) {
      long blockShift = (i*blocksPerCtxt+j)*16;
      for (long k=0; k<16; k++) {
        long byteIdx = blockShift+k;
        if (byteIdx < data.length()) {
          long slotIdx = j + k*blocksPerCtxt;
          GF2XFromBytes(slots[slotIdx], &data[byteIdx], 1);
        }
      }
    }
    ea2.encode(encData[i], slots);
  }
}

void decode4AES(Vec<uint8_t>& data, const Vec<ZZX>& encData,
                const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long nBytes = encData.length() * ea2.size();
  if (data.length()<=0 || data.length()>nBytes) data.SetLength(nBytes);
  long nBlocks = helib::divc(data.length(), 16);
  long blocksPerCtxt = ea2.size() / 16;
  long nCtxt = helib::divc(nBlocks, blocksPerCtxt);
  vector<GF2X> slots;
  for (long i=0; i<nCtxt; i++) {
    ea2.decode(slots, encData[i]);
    for (long j=0; j<blocksPerCtxt; j++) {
      long blockShift = (i*blocksPerCtxt+j)*16;
      for (long k=0; k<16; k++) {
        long byteIdx = blockShift+k;
        if (byteIdx < data.length()) {
          long slotIdx = j + k*blocksPerCtxt;
          BytesFromGF2X(&data[byteIdx], slots[slotIdx], 1);
        }
      }
    }
  }
}

static void invertOne(helib::Ctxt& c, bool verbose)
{
  if(verbose) std::cerr << "  [invert] start noise log2=" << log(c.getNoiseBound())/log(2.0) << std::endl;
  helib::Ctxt tmp1(c);
  tmp1.frobeniusAutomorph(1);
  if(verbose) std::cerr << "  [invert] after frob(1) on tmp1: bitCap=" << tmp1.bitCapacity()
                      << " noise log2=" << log(tmp1.getNoiseBound())/log(2.0) << std::endl;
  c.multiplyBy(tmp1);
  if(verbose) std::cerr << "  [invert] after mul1 (data*=tmp1): bitCap=" << c.bitCapacity()
                      << " noise log2=" << log(c.getNoiseBound())/log(2.0) << std::endl;
  helib::Ctxt tmp2(c);
  tmp2.frobeniusAutomorph(2);
  tmp1.multiplyBy(tmp2);
  if(verbose) std::cerr << "  [invert] after mul2: data bitCap=" << c.bitCapacity()
                      << " tmp1 bitCap=" << tmp1.bitCapacity()
                      << " data noise log2=" << log(c.getNoiseBound())/log(2.0)
                      << " tmp1 noise log2=" << log(tmp1.getNoiseBound())/log(2.0) << std::endl;
  c.multiplyBy(tmp2);
  if(verbose) std::cerr << "  [invert] after mul3: bitCap=" << c.bitCapacity()
                      << " noise log2=" << log(c.getNoiseBound())/log(2.0) << std::endl;
  c.frobeniusAutomorph(4);
  c.multiplyBy(tmp1);
  if(verbose) std::cerr << "  [invert] after mul4: bitCap=" << c.bitCapacity()
                      << " noise log2=" << log(c.getNoiseBound())/log(2.0) << std::endl;
}
static void invert(vector<helib::Ctxt>& data)
{
  const size_t maxThreads = std::min<size_t>(64, data.size());
  std::vector<std::thread> threads;
  std::atomic<size_t> nextIdx(0);

  auto worker = [&]() {
    size_t i;
    while ((i = nextIdx.fetch_add(1)) < data.size()) {
      invertOne(data[i], (i==0));
    }
  };

  for (size_t t = 0; t < maxThreads; t++) {
    threads.emplace_back(worker);
  }
  for (auto& t : threads) t.join();
}

static void packCtxt(vector<helib::Ctxt>& to, const vector<helib::Ctxt>& from,
                     const GF2X& XinSlots)
{
  HELIB_TIMER_START;
  if (from.size() <= 1) { to = from; return; }
  const helib::Context& context = from[0].getContext();
  const helib::EncryptedArrayDerived<helib::PA_GF2>& ea =
    context.getEA().getDerived(helib::PA_GF2());
  const GF2XModulus& PhimX = ea.getTab().getPhimXMod();
  long e = ea.getDegree() / 8;
  long nPacked = helib::divc(from.size(), e);
  to.assign(nPacked, helib::Ctxt(helib::ZeroCtxtLike, from[0]));
  auto _pkBody = [&](long i) {
    to[i] = from[e*i];
    if (e*i == (long)from.size()-1) return;
    helib::Ctxt tmp = from[e*i+1];
    tmp.multByConstant(conv<ZZX>(XinSlots)); to[i] += tmp;
    if (e*i == (long)from.size()-2) return;
    GF2X X2i = XinSlots;
    for (long j=e*i+2; j<e*(i+1) && j<(long)from.size(); j++) {
      MulMod(X2i, X2i, XinSlots, PhimX);
      tmp = from[j]; tmp.multByConstant(conv<ZZX>(X2i)); to[i] += tmp;
    }
  };
  if ((long)to.size() < 8) {
    for (long i=0; i<(long)to.size(); i++) _pkBody(i);
  } else {
    const size_t _pkMax = std::min<size_t>(64, (size_t)to.size());
    std::vector<std::thread> _pkThreads;
    std::atomic<size_t> _pkNext(0);
    auto _pkWorker = [&]() {
      size_t idxS;
      while ((idxS = _pkNext.fetch_add(1)) < (size_t)to.size()) _pkBody((long)idxS);
    };
    for (size_t t = 0; t < _pkMax; t++) _pkThreads.emplace_back(_pkWorker);
    for (auto& t : _pkThreads) t.join();
  }
}

namespace helib {
class AESHoistPrecon
{
  helib::Ctxt ctxt;
  NTL::xdouble noise;
  std::vector<helib::DoubleCRT> polyDigits;

public:
  AESHoistPrecon(const helib::Ctxt& _ctxt) : ctxt(_ctxt), noise(1.0)
  {
    if (ctxt.parts.size() >= 1)
      helib::assertTrue(
          ctxt.parts[0].skHandle.isOne(),
          "Invalid ciphertext (secret key handle for part 0 is not one)");
    if (ctxt.parts.size() <= 1)
      return;

    ctxt.cleanUp();
    const helib::Context& context = ctxt.getContext();
    const helib::PubKey& pubKey = ctxt.getPubKey();
    long keyID = ctxt.getKeyID();

    helib::assertTrue(ctxt.inCanonicalForm(keyID),
               "Ciphertext is not in canonical form");

    ctxt.relin_CKKS_adjust();

    NTL::xdouble addedNoise = ctxt.parts[1].breakIntoDigits(polyDigits);
    NTL::xdouble max_ks_noise(0.0);
    for (const helib::KeySwitch& ks : pubKey.keySWlist()) {
      if (max_ks_noise < ks.noiseBound)
        max_ks_noise = ks.noiseBound;
    }
    addedNoise *= max_ks_noise;

    double logProd = context.logOfProduct(context.getSpecialPrimes());
    noise = ctxt.getNoiseBound() * NTL::xexp(logProd);

    double ratio = NTL::conv<double>(addedNoise / noise);
    if (ratio > 1) {
      helib::Warning("KS-noise-ratio-hoist=" + std::to_string(ratio));
    }
    noise += addedNoise;
  }

  std::shared_ptr<helib::Ctxt> automorph(long k) const
  {
    std::shared_ptr<helib::Ctxt> result = std::make_shared<helib::Ctxt>(ctxt);

    if (k == 1 || ctxt.isEmpty())
      return result;

    const helib::Context& context = ctxt.getContext();
    const helib::PubKey& pubKey = ctxt.getPubKey();

    result = std::make_shared<helib::Ctxt>(helib::ZeroCtxtLike, ctxt);
    result->noiseBound = noise;
    result->intFactor = ctxt.intFactor;
    result->primeSet = ctxt.primeSet | context.getSpecialPrimes();

    if (ctxt.parts.size() == 1) {
      helib::CtxtPart tmpPart = ctxt.parts[0];
      tmpPart.automorph(k);
      tmpPart.addPrimesAndScale(context.getSpecialPrimes());
      result->addPart(tmpPart, true);
      return result;
    }

    long keyID = ctxt.getKeyID();
    if (!pubKey.isReachable(k, keyID)) {
      throw helib::LogicError("no key-switching matrices for k=" + std::to_string(k) +
                       ", keyID=" + std::to_string(keyID));
    }

    const helib::KeySwitch& W = pubKey.getNextKSWmatrix(k, keyID);
    long amt = W.fromKey.getPowerOfX();

    helib::CtxtPart tmpPart = ctxt.parts[0];
    tmpPart.automorph(amt);
    tmpPart.addPrimesAndScale(context.getSpecialPrimes());
    result->addPart(tmpPart, true);

    std::vector<helib::DoubleCRT> tmpDigits = polyDigits;
    for (auto&& tmp : tmpDigits)
      tmp.automorph(amt);
    result->keySwitchDigits(W, tmpDigits);
    long m = context.getM();
    if ((amt - k) % m != 0) {
      long k2 = NTL::MulMod(k, NTL::InvMod(amt, m), m);
      result->smartAutomorph(k2);
    }
    return result;
  }
};
} // namespace helib
static void unackCtxt(vector<helib::Ctxt>& to, const vector<helib::Ctxt>& from,
                      const Mat<GF2X>& unpackConsts)
{
  HELIB_TIMER_START;
  const helib::Context& context = from[0].getContext();
  const helib::EncryptedArrayDerived<helib::PA_GF2>& ea =
    context.getEA().getDerived(helib::PA_GF2());
  long e = ea.getDegree() / 8;
  long nUnpacked = from.size()*e;
  if (to.size()==0) to.resize(nUnpacked, helib::Ctxt(helib::ZeroCtxtLike, from[0]));
  else {
    if (nUnpacked > (long)to.size()) nUnpacked = to.size();
    for (long i=0; i<(long)to.size(); i++) to[i].clear();
  }
  long nPacked = helib::divc(nUnpacked, e);
  auto _unBody = [&](long idx) {
    vector<helib::Ctxt> conjugates(e, from[idx]);
    for (long j=1; j<e; j++) conjugates[j].frobeniusAutomorph(8*j);
    for (long i=0; i<e && (idx*e+i)<(long)to.size(); i++) {
      for (long j=0; j<e; j++) {
        helib::Ctxt tmp = conjugates[j];
        tmp.multByConstant(conv<ZZX>(unpackConsts[i][j]));
        to[idx*e+i] += tmp;
      }
    }
  };
  if (nPacked < 8) {
    for (long idx=0; idx<nPacked; idx++) _unBody(idx);
  } else {
    const size_t _unMax = std::min<size_t>(64, (size_t)nPacked);
    std::vector<std::thread> _unThreads;
    std::atomic<size_t> _unNext(0);
    auto _unWorker = [&]() {
      size_t idxS;
      while ((idxS = _unNext.fetch_add(1)) < (size_t)nPacked) _unBody((long)idxS);
    };
    for (size_t t = 0; t < _unMax; t++) _unThreads.emplace_back(_unWorker);
    for (auto& t : _unThreads) t.join();
  }
}

// =========================================================================
// 方案A 第0步验证: byte -> 8bit -> byte 往返链路
// =========================================================================
static void buildByteToBitMat(vector<PolyType>& bitMat,
                               const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  bitMat.clear();
  bitMat.resize(8*8, helib::DoubleCRT(0, ea2.getContext(),
                helib::IndexSet(0, ea2.getContext().numPrimes()-1)));

  vector<GF2X> encScratch(ea2.size());
  GF2X zero, one; SetCoeff(one, 0, 1);

  for (long k=0; k<8; k++) {
    vector<GF2X> mapOut(8);
    for (long i=0;i<8;i++) mapOut[i] = (i==k) ? one : zero;

    vector<GF2X> Ck;
    ea2.buildLinPolyCoeffs(Ck, mapOut);

    for (long j=0; j<8; j++) {
      for (long s=0; s<ea2.size(); s++) encScratch[s] = Ck[j];
      ZZX tmpZZX; ea2.encode(tmpZZX, encScratch);
      bitMat[k*8+j] = tmpZZX;
    }
  }
}


static void unpackByteToBits(vector<helib::Ctxt>& bits, const helib::Ctxt& byteCtxt,
                              const vector<PolyType>& bitMat,
                              const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  long d = ea2.getDegree();

  vector<helib::Ctxt> conjugates(d, byteCtxt);
  {
    const helib::Context& hctx = byteCtxt.getContext();
    long m_val = hctx.getM();
    long p_val = hctx.getP();
    long d_val = hctx.getOrdP();
    helib::AESHoistPrecon precon(byteCtxt);
    for (long j=1; j<d; j++) {
      long jj = helib::mcMod(j, d_val);
      long val = NTL::PowerMod(p_val % m_val, jj, m_val);
      auto res = precon.automorph(val);
      conjugates[j] = *res;
    }
  }

  bits.assign(8, byteCtxt);
  for (long k=0; k<8; k++) {
    helib::Ctxt acc = conjugates[0];
    acc.multByConstant(bitMat[k*8+0]);
    for (long j=1; j<d; j++) {
      helib::Ctxt term = conjugates[j];
      term.multByConstant(bitMat[k*8+j]);
      acc += term;
    }
    bits[k] = acc;
  }
}

static void packBitsToByte(helib::Ctxt& byteCtxt, const vector<helib::Ctxt>& bits,
                            const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2)
{
  vector<GF2X> encScratch(ea2.size());
  byteCtxt = bits[0];

  for (long k=1; k<8; k++) {
    helib::Ctxt term = bits[k];
    GF2X zetaK(k,1);
    for (long s=0; s<ea2.size(); s++) encScratch[s] = zetaK;
    ZZX constPoly; ea2.encode(constPoly, encScratch);
    term.multByConstant(constPoly);
    byteCtxt += term;
  }
}


// combinedLinearLayer: 把 affine(bit层面) + ShiftRows + MixColumns 合并成一次128-bit线性变换
// 输入/输出: stateBits[j][k], j=0..15 字节位置(idx(r,c)=r+4c), k=0..7 bit位置(bit0=LSB)
static void combinedLinearLayer(vector<vector<helib::Ctxt>>& stateBits,
                                 const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2,
                                 bool isLastRound)
{
  // stateBits.size() 现在可以是 numBlocks*16 (每个block独立16字节=128bit,
  // 各自应用同一套128x128合并线性变换, block之间互不干扰)
  long numBlocks = (long)stateBits.size() / 16;
  helib::assertTrue<helib::InvalidArgument>((long)stateBits.size() == numBlocks*16,
    "stateBits.size() must be a multiple of 16");
  vector<helib::Ctxt> flat;
  flat.reserve((size_t)numBlocks*128);
  for (long b=0; b<numBlocks; b++)
    for (long j=0; j<16; j++)
      for (long k=0; k<8; k++)
        flat.push_back(stateBits[b*16+j][k]);
  vector<GF2X> encScratch(ea2.size());
  GF2X one; SetCoeff(one, 0, 1);
  for (long s=0; s<ea2.size(); s++) encScratch[s] = one;
  ZZX constOnePoly; ea2.encode(constOnePoly, encScratch);
  size_t totalOut = (size_t)numBlocks*128;
  vector<helib::Ctxt> out(totalOut, flat[0]);
  {
    const size_t _lMax = std::min<size_t>(32, totalOut);
    std::vector<std::thread> _lThreads;
    std::atomic<size_t> _lNext(0);
    auto _lWorker = [&]() {
      size_t idx;
      while ((idx = _lNext.fetch_add(1)) < totalOut) {
        long b = (long)(idx / 128);
        int i = (int)(idx % 128);
        long base = b*128;
        const std::vector<int>& terms = isLastRound ? COMBINED_LIN_ROWS_LASTROUND[i] : COMBINED_LIN_ROWS[i];
        int constBit = isLastRound ? COMBINED_LIN_CONST_LASTROUND[i] : COMBINED_LIN_CONST[i];
        helib::Ctxt acc = flat[base+terms[0]];
        for (size_t t=1; t<terms.size(); t++) acc += flat[base+terms[t]];
        if (constBit) acc.addConstant(constOnePoly);
        out[idx] = acc;
      }
    };
    for (size_t t = 0; t < _lMax; t++) _lThreads.emplace_back(_lWorker);
    for (auto& t : _lThreads) t.join();
  }
  for (long b=0; b<numBlocks; b++)
    for (long j=0; j<16; j++)
      for (long k=0; k<8; k++)
        stateBits[b*16+j][k] = out[(size_t)b*128 + j*8+k];
}


// combinedLinearLayerDecrypt: InvMixColumns + InvShiftRows + InvSubBytes仿射部分 合并成一次128-bit线性变换
// isFirstRound: 第一轮(紧接最初AddRoundKey之后)不含InvMixColumns
static void combinedLinearLayerDecrypt(vector<vector<helib::Ctxt>>& stateBits,
                                        const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2,
                                        bool isFirstRound)
{
  vector<helib::Ctxt> flat;
  flat.reserve(128);
  for (long j=0; j<16; j++)
    for (long k=0; k<8; k++)
      flat.push_back(stateBits[j][k]);

  vector<GF2X> encScratch(ea2.size());
  GF2X one; SetCoeff(one, 0, 1);
  for (long s=0; s<ea2.size(); s++) encScratch[s] = one;
  ZZX constOnePoly; ea2.encode(constOnePoly, encScratch);

  vector<helib::Ctxt> out;
  out.reserve(128);
  for (int i=0; i<128; i++) {
    const std::vector<int>& terms = isFirstRound ? DEC_COMBINED_LIN_ROWS_FIRSTROUND[i] : DEC_COMBINED_LIN_ROWS[i];
    int constBit = isFirstRound ? DEC_COMBINED_LIN_CONST_FIRSTROUND[i] : DEC_COMBINED_LIN_CONST[i];
    helib::Ctxt acc = flat[terms[0]];
    for (size_t t=1; t<terms.size(); t++) acc += flat[terms[t]];
    if (constBit) acc.addConstant(constOnePoly);
    out.push_back(acc);
  }

  for (long j=0; j<16; j++)
    for (long k=0; k<8; k++)
      stateBits[j][k] = out[j*8+k];
}

void HomAES::debugUnpackTest(const helib::PubKey& pk, const helib::SecKey& sk) const
{
  cout << "\\n=== [方案A验证] unpack/pack byte<->bit 链路测试 ===" << endl;

  vector<PolyType> bitMat;
  buildByteToBitMat(bitMat, ea2);
  cout << "buildByteToBitMat 完成" << endl;

  // 构造一个已知测试byte
  uint8_t testByte = 0xA5; // 10100101
  cout << "测试byte = 0x" << hex << (int)testByte << dec << endl;
  cout << "真实8个bit(从bit0/LSB到bit7/MSB): ";
  for (int b=0;b<8;b++) cout << ((testByte>>b)&1) << " ";
  cout << endl;

  // 构造一个整个ctxt都填同一个byte的明文向量(简单起见,全部960 slot都填testByte对应的域元素)
  Vec<uint8_t> testVec(INIT_SIZE, ea2.size()/16*16);
  for (long i=0; i<testVec.length(); i++) testVec[i] = testByte;

  Vec<ZZX> encoded;
  encode4AES(encoded, testVec, ea2);
  cout << "encode4AES 完成, encoded.length()=" << encoded.length() << endl;

  // 加密成一个Ctxt (只用第一个, 因为16个byte位置的encoded都一样内容,取第0个即可)
  helib::Ctxt byteCtxt(pk);
  pk.Encrypt(byteCtxt, encoded[0]);
  cout << "加密完成, bitCap=" << byteCtxt.bitCapacity() << endl;

  // unpack
  vector<helib::Ctxt> bits;
  unpackByteToBits(bits, byteCtxt, bitMat, ea2);
  cout << "unpackByteToBits 完成, bits.size()=" << bits.size()
       << " bitCap[0]=" << bits[0].bitCapacity() << endl;

  // 解密每个bit密文, 检查slot0(对应第一个AES block第一个byte位置)是否等于真实bit
  cout << "解密验证8个bit密文: ";
  bool unpackOK = true;
  for (int k=0; k<8; k++) {
    ZZX poly;
    sk.Decrypt(poly, bits[k]);
    vector<GF2X> slots;
    ea2.decode(slots, poly);
    // slots[0] 应该是常数0或常数1 (GF2X表示的0次多项式)
    long bitVal = IsZero(slots[0]) ? 0 : 1;
    long expected = (testByte >> k) & 1;
    cout << bitVal << "(期望" << expected << ") ";
    if (bitVal != expected) unpackOK = false;
  }
  cout << endl;
  cout << (unpackOK ? "unpack验证通过 ✓" : "unpack验证失败 ✗") << endl;

  // pack回byte
  helib::Ctxt repacked(pk);
  packBitsToByte(repacked, bits, ea2);
  cout << "packBitsToByte 完成, bitCap=" << repacked.bitCapacity() << endl;

  ZZX repackedPoly;
  sk.Decrypt(repackedPoly, repacked);
  vector<GF2X> repackedSlots;
  ea2.decode(repackedSlots, repackedPoly);

  // 把repackedSlots[0]转回byte比较 (用跟encode4AES对应的方式还原, 简化：直接转成字节数组比较低8位)
  Vec<uint8_t> decodedBack(INIT_SIZE, 1);
  BytesFromGF2X(decodedBack.data(), repackedSlots[0], 1);
  cout << "pack回去的byte = 0x" << hex << (int)decodedBack[0] << dec
       << " (期望 0x" << hex << (int)testByte << dec << ")" << endl;

  bool packOK = (decodedBack[0] == testByte);
  cout << (packOK ? "pack往返验证通过 ✓" : "pack往返验证失败 ✗") << endl;

  cout << "=== 验证" << ((unpackOK && packOK) ? "全部通过 ✓✓✓" : "存在失败 ✗") << " ===" << endl;
}

void HomAES::batchRecryptIsolationTest(const helib::PubKey& pk, const helib::SecKey& sk) const
{
  cout << "\n=== [隔离测试] batchRecrypt 是否保值 ===" << endl;

  // 用一个已知的、非平凡的测试值(不用全0全1, 避免掩盖bug)
  uint8_t testBlock[16];
  for (int i = 0; i < 16; i++) testBlock[i] = (uint8_t)(0x10 * i + 0x03);  // 03 13 23 33 ... f3

  cout << "[隔离测试] 原始测试值: ";
  for (int i = 0; i < 16; i++) fprintf(stdout, "%02x ", testBlock[i]);
  cout << endl;

  vector<ZZX> encoded16;
  encode4AES_singleBlock16Ctxt(encoded16, testBlock, ea2);

  vector<helib::Ctxt> eDataTest;
  eDataTest.reserve(16);
  for (long j = 0; j < 16; j++) {
    eDataTest.push_back(helib::Ctxt(pk));
    pk.Encrypt(eDataTest[j], encoded16[j]);
  }

  cout << "[隔离测试] 加密完成, bitCap=" << eDataTest[0].bitCapacity() << endl;

  // 解密一次, 确认加密本身没问题 (作为对照组基线)
  {
    vector<ZZX> polyCheck(16);
    for (long i = 0; i < 16; i++) sk.Decrypt(polyCheck[i], eDataTest[i]);
    uint8_t checkBlock[16];
    decode4AES_singleBlock16Ctxt(checkBlock, polyCheck, ea2);
    bool matchBefore = true;
    for (int i = 0; i < 16; i++) if (checkBlock[i] != testBlock[i]) matchBefore = false;
    cout << "[隔离测试] batchRecrypt之前解密验证: " << (matchBefore ? "一致 ✓" : "不一致 ✗ (加密/编码本身就有问题)") << endl;
  }

  // 核心测试: 调用 batchRecrypt
  double t0 = NTL::GetTime();
  batchRecrypt(eDataTest);
  double t1 = NTL::GetTime();
  cout << "[隔离测试] batchRecrypt 耗时=" << (t1 - t0) << "s, 之后 bitCap=" << eDataTest[0].bitCapacity() << endl;

  // 解密, 检查 batchRecrypt 之后数值是否还保持一致
  vector<ZZX> polyAfter(16);
  for (long i = 0; i < 16; i++) sk.Decrypt(polyAfter[i], eDataTest[i]);
  uint8_t afterBlock[16];
  decode4AES_singleBlock16Ctxt(afterBlock, polyAfter, ea2);

  cout << "[隔离测试] batchRecrypt之后的值: ";
  for (int i = 0; i < 16; i++) fprintf(stdout, "%02x ", afterBlock[i]);
  cout << endl;

  bool matchAfter = true;
  for (int i = 0; i < 16; i++) if (afterBlock[i] != testBlock[i]) matchAfter = false;

  if (matchAfter) {
    cout << "[隔离测试] batchRecrypt 保值 ✓ ✓ ✓" << endl;
  } else {
    cerr << "[隔离测试] batchRecrypt 破坏了数据 ✗ ✗ ✗ (确认是 packCtxt/unackCtxt 打包格式不兼容问题)" << endl;
  }
}

// =========================================================================
// 方案A 完整10轮: invert(复用) -> unpack -> bit级(ShiftRows+MixColumns) -> pack -> +key
// =========================================================================

static void bitLevelShiftRows(vector<vector<helib::Ctxt>>& stateBits)
{
  vector<vector<helib::Ctxt>> tmp = stateBits;
  auto idx = [](int r,int c){ return r + 4*c; };
  for (int c=0;c<4;c++) stateBits[idx(0,c)] = tmp[idx(0,c)];
  for (int c=0;c<4;c++) stateBits[idx(1,c)] = tmp[idx(1,(c+1)%4)];
  for (int c=0;c<4;c++) stateBits[idx(2,c)] = tmp[idx(2,(c+2)%4)];
  for (int c=0;c<4;c++) stateBits[idx(3,c)] = tmp[idx(3,(c+3)%4)];
}

static const int M2_bits[8][8] = {
    {0,0,0,0,0,0,0,1},{1,0,0,0,0,0,0,1},{0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,1},
    {0,0,0,1,0,0,0,1},{0,0,0,0,1,0,0,0},{0,0,0,0,0,1,0,0},{0,0,0,0,0,0,1,0}
};
static void matVecMulCtxt(const int M[8][8], const vector<helib::Ctxt>& in,
                           vector<helib::Ctxt>& out)
{
  if (in.size() != 8) {
    std::cerr << "[FATAL] matVecMulCtxt: in.size()=" << in.size()
               << " (expected 8)" << std::endl;
    std::abort();
  }
  out.clear(); out.reserve(8);
  for (int i=0;i<8;i++) {
    bool first = true;
    helib::Ctxt acc = in[0];
    for (int j=0;j<8;j++) {
      if (M[i][j]) {
        if (first) { acc = in[j]; first=false; }
        else acc += in[j];
      }
    }
    if (first) acc.clear();
    out.push_back(acc);
  }
}
static void bitLevelMixColumns(vector<vector<helib::Ctxt>>& stateBits)
{
  static int M3_bits[8][8]; static bool inited=false;
  if(!inited){ for(int i=0;i<8;i++) for(int j=0;j<8;j++) M3_bits[i][j]=M2_bits[i][j]^(i==j?1:0); inited=true; }
  auto idx = [](int r,int c){ return r + 4*c; };
  for (int c=0;c<4;c++) {
    vector<helib::Ctxt> a0=stateBits[idx(0,c)], a1=stateBits[idx(1,c)],
                         a2=stateBits[idx(2,c)], a3=stateBits[idx(3,c)];
    vector<helib::Ctxt> o0,o1,o2,o3,tmp;
    matVecMulCtxt(M2_bits, a0, o0);
    matVecMulCtxt(M3_bits, a1, tmp); for(int i=0;i<8;i++) o0[i]+=tmp[i];
    for(int i=0;i<8;i++) o0[i]+=a2[i];
    for(int i=0;i<8;i++) o0[i]+=a3[i];

    o1 = a0;
    matVecMulCtxt(M2_bits, a1, tmp); for(int i=0;i<8;i++) o1[i]+=tmp[i];
    matVecMulCtxt(M3_bits, a2, tmp); for(int i=0;i<8;i++) o1[i]+=tmp[i];
    for(int i=0;i<8;i++) o1[i]+=a3[i];

    o2 = a0; for(int i=0;i<8;i++) o2[i]+=a1[i];
    matVecMulCtxt(M2_bits, a2, tmp); for(int i=0;i<8;i++) o2[i]+=tmp[i];
    matVecMulCtxt(M3_bits, a3, tmp); for(int i=0;i<8;i++) o2[i]+=tmp[i];

    matVecMulCtxt(M3_bits, a0, o3);
    for(int i=0;i<8;i++) o3[i]+=a1[i];
    for(int i=0;i<8;i++) o3[i]+=a2[i];
    matVecMulCtxt(M2_bits, a3, tmp); for(int i=0;i<8;i++) o3[i]+=tmp[i];

    stateBits[idx(0,c)]=o0; stateBits[idx(1,c)]=o1;
    stateBits[idx(2,c)]=o2; stateBits[idx(3,c)]=o3;
  }
}

void HomAES::testSboxOnly(vector<helib::Ctxt>& eData) const
{
  invert(eData);
  for (long j=0; j<(long)eData.size(); j++) {
    helib::applyLinPolyLL(eData[j], encAffMat, ea2.getDegree());
    eData[j].addConstant(affVec);
  }
}
void HomAES::bitLevelRoundTestA(vector<helib::Ctxt>& eData,
                                 const vector<vector<helib::Ctxt>>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return;
  static vector<PolyType> bitMat;
  if (bitMat.empty()) buildByteToBitMat(bitMat, ea2);
  for (long k=0; k<(long)eData.size(); k++) eData[k] += aesKey[0][k];
  double t_invert=0, t_unpack=0, t_linear=0, t_pack=0, t_boot=0;
  double wall_boot=0, wall_invert=0, wall_unpack=0, wall_linear=0, wall_pack=0;
  std::chrono::steady_clock::time_point _wnow0;
  long bootCount = 0;
  for (long i=1; i<(long)aesKey.size(); i++) {
    double t0, t1;
    bool isLastRound = (i == (long)aesKey.size()-1);
    std::cerr << "[A pre-invert] bitCap=" << eData[0].bitCapacity() << " round=" << i << std::endl;
    if (eData[0].bitCapacity() < 80) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] ENC round=" << i << " (pre-invert, thresh80) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
    }
    std::cerr << "  [noise] before invert: " << eData[0].getNoiseBound()
               << " log2=" << log(eData[0].getNoiseBound())/log(2.0) << std::endl;
    auto _wall0 = std::chrono::steady_clock::now();
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); invert(eData); t1=NTL::GetTime(); t_invert+=t1-t0;
    auto _wall1 = std::chrono::steady_clock::now();
    double _wallSec = std::chrono::duration<double>(_wall1-_wall0).count();
    wall_invert += _wallSec;
    std::cerr << "  [WALL] invert_wall=" << _wallSec << "s (NTL_time=" << (t1-t0) << "s)" << std::endl;
    std::cerr << "  [A TIME] round=" << i << " invert=" << t1-t0
               << "s bitCap=" << eData[0].bitCapacity() << std::endl;
    std::cerr << "  [noise] after invert: " << eData[0].getNoiseBound()
               << " log2=" << log(eData[0].getNoiseBound())/log(2.0) << std::endl;
    if (eData[0].bitCapacity() < 35) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] ENC round=" << i << " (pre-unpack, thresh35) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
    }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    vector<vector<helib::Ctxt>> stateBits(eData.size());
    {
      std::atomic<bool> _unpackFatal(false);
      std::atomic<long> _unpackFatalIdx(-1);
      const size_t _uMax = std::min<size_t>(64, eData.size());
      std::vector<std::thread> _uThreads;
      std::atomic<size_t> _uNext(0);
      auto _uWorker = [&]() {
        size_t j;
        while ((j = _uNext.fetch_add(1)) < eData.size()) {
          unpackByteToBits(stateBits[j], eData[j], bitMat, ea2);
          if (stateBits[j].size() != 8) {
            _unpackFatal = true;
            _unpackFatalIdx = (long)j;
          }
        }
      };
      for (size_t t = 0; t < _uMax; t++) _uThreads.emplace_back(_uWorker);
      for (auto& t : _uThreads) t.join();
      if (_unpackFatal) {
        std::cerr << "[FATAL] unpackByteToBits produced stateBits[" << _unpackFatalIdx
                   << "].size() != 8 (expected 8)" << std::endl;
        std::abort();
      }
    }
    t1=NTL::GetTime(); t_unpack += t1-t0; wall_unpack+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A TIME] round=" << i << " unpack=" << t1-t0
               << "s bitCap=" << stateBits[0][0].bitCapacity() << std::endl;
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    combinedLinearLayer(stateBits, ea2, isLastRound);
    t1=NTL::GetTime(); t_linear += t1-t0; wall_linear+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A TIME] round=" << i << " combinedLinear(affine+shift+mix)=" << t1-t0
               << "s bitCap=" << stateBits[0][0].bitCapacity() << std::endl;
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    {
      const size_t _pMax = std::min<size_t>(64, eData.size());
      std::vector<std::thread> _pThreads;
      std::atomic<size_t> _pNext(0);
      auto _pWorker = [&]() {
        size_t j;
        while ((j = _pNext.fetch_add(1)) < eData.size()) {
          packBitsToByte(eData[j], stateBits[j], ea2);
        }
      };
      for (size_t t = 0; t < _pMax; t++) _pThreads.emplace_back(_pWorker);
      for (auto& t : _pThreads) t.join();
    }
    t1=NTL::GetTime(); t_pack += t1-t0; wall_pack+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A TIME] round=" << i << " pack=" << t1-t0
               << "s bitCap=" << eData[0].bitCapacity() << std::endl;
    for (long k=0; k<(long)eData.size(); k++) eData[k] += aesKey[i][k];
    std::cerr << "  [A TIME] round=" << i << " after addkey bitCap="
               << eData[0].bitCapacity() << std::endl;
  }
  std::cerr << "[A TOTAL] invert=" << t_invert << "s unpack=" << t_unpack
             << "s combinedLinear=" << t_linear << "s pack=" << t_pack
             << "s boot=" << t_boot << "s bootCount=" << bootCount << std::endl;
  std::cerr << "[A TOTAL][WALL] invert_wall=" << wall_invert
             << "s unpack_wall=" << wall_unpack << "s linear_wall=" << wall_linear
             << "s pack_wall=" << wall_pack << "s boot_wall=" << wall_boot << "s" << std::endl;
  double total_ntl = t_invert + t_unpack + t_linear + t_pack + t_boot;
  double total_wall = wall_invert + wall_unpack + wall_linear + wall_pack + wall_boot;
  std::cerr << "[A TOTAL] total=" << total_ntl << "s" << std::endl;
  std::cerr << "[A TOTAL][WALL] total_wall=" << total_wall << "s" << std::endl;
}

void HomAES::bitLevelRoundTestADecrypt(vector<helib::Ctxt>& eData,
                                        const vector<vector<helib::Ctxt>>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return;
  static vector<PolyType> bitMat;
  if (bitMat.empty()) buildByteToBitMat(bitMat, ea2);
  long nRounds = (long)aesKey.size();  // 11 (含初始轮密钥)
  for (long k=0; k<(long)eData.size(); k++) eData[k] += aesKey[nRounds-1][k];
  double t_invert=0, t_unpack=0, t_linear=0, t_pack=0, t_boot=0;
  double wall_boot=0, wall_invert=0, wall_unpack=0, wall_linear=0, wall_pack=0;
  std::chrono::steady_clock::time_point _wnow0;
  long bootCount = 0;
  for (long rnd=nRounds-2; rnd>=0; rnd--) {
    double t0, t1;
    bool isFirstIter = (rnd == nRounds-2);

    std::cerr << "[A-DEC pre-unpack] bitCap=" << eData[0].bitCapacity() << " rnd=" << rnd << std::endl;
    if (eData[0].bitCapacity() < 35) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] DEC rnd=" << rnd << " (pre-unpack, thresh35) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
    }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    vector<vector<helib::Ctxt>> stateBits(16);
    for (long j=0; j<(long)eData.size(); j++) {
      unpackByteToBits(stateBits[j], eData[j], bitMat, ea2);
      if (stateBits[j].size() != 8) {
        std::cerr << "[FATAL-DEC] unpackByteToBits produced stateBits[" << j
                   << "].size()=" << stateBits[j].size() << " (expected 8)" << std::endl;
        std::abort();
      }
    }
    t1=NTL::GetTime(); t_unpack += t1-t0; wall_unpack+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A-DEC TIME] rnd=" << rnd << " unpack=" << t1-t0
               << "s bitCap=" << stateBits[0][0].bitCapacity() << std::endl;

    std::cerr << "[A-DEC pre-linear] bitCap=" << stateBits[0][0].bitCapacity() << " rnd=" << rnd << std::endl;
    if (stateBits[0][0].bitCapacity() < 35) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
      for (long j=0; j<(long)eData.size(); j++)
        packBitsToByte(eData[j], stateBits[j], ea2);
      auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] DEC rnd=" << rnd << " (pre-linear, thresh35, incl pack) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
      for (long j=0; j<(long)eData.size(); j++)
        unpackByteToBits(stateBits[j], eData[j], bitMat, ea2);
    }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    combinedLinearLayerDecrypt(stateBits, ea2, isFirstIter);
    t1=NTL::GetTime(); t_linear += t1-t0; wall_linear+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A-DEC TIME] rnd=" << rnd << " combinedLinearDecrypt(invmix+invshift+invaffine)=" << t1-t0
               << "s bitCap=" << stateBits[0][0].bitCapacity() << std::endl;

    std::cerr << "[A-DEC pre-pack] bitCap=" << stateBits[0][0].bitCapacity() << " rnd=" << rnd << std::endl;
    if (stateBits[0][0].bitCapacity() < 35) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
      for (long j=0; j<(long)eData.size(); j++)
        packBitsToByte(eData[j], stateBits[j], ea2);
      auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] DEC rnd=" << rnd << " (pre-pack, thresh35, incl pack) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
      for (long j=0; j<(long)eData.size(); j++)
        unpackByteToBits(stateBits[j], eData[j], bitMat, ea2);
    }
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now();
    for (long j=0; j<(long)eData.size(); j++)
      packBitsToByte(eData[j], stateBits[j], ea2);
    t1=NTL::GetTime(); t_pack += t1-t0; wall_pack+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_wnow0).count();
    std::cerr << "  [A-DEC TIME] rnd=" << rnd << " pack=" << t1-t0
               << "s bitCap=" << eData[0].bitCapacity() << std::endl;

    std::cerr << "[A-DEC pre-invert] bitCap=" << eData[0].bitCapacity() << " rnd=" << rnd << std::endl;
    if (eData[0].bitCapacity() < 100) {
      t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); auto _bw0=std::chrono::steady_clock::now(); batchRecrypt(eData); double _bt=NTL::GetTime()-t0; t_boot+=_bt; bootCount++; wall_boot+=std::chrono::duration<double>(std::chrono::steady_clock::now()-_bw0).count();
      std::cerr << "  [BOOT] DEC rnd=" << rnd << " (pre-invert, thresh100) elapsed=" << _bt
                 << "s bitCap_after=" << eData[0].bitCapacity() << std::endl;
    }
    std::cerr << "  [noise] before invert(dec): " << eData[0].getNoiseBound()
               << " log2=" << log(eData[0].getNoiseBound())/log(2.0) << std::endl;
    auto _wall0 = std::chrono::steady_clock::now();
    t0=NTL::GetTime(); _wnow0=std::chrono::steady_clock::now(); invert(eData); t1=NTL::GetTime(); t_invert+=t1-t0;
    auto _wall1 = std::chrono::steady_clock::now();
    double _wallSec = std::chrono::duration<double>(_wall1-_wall0).count();
    wall_invert += _wallSec;
    std::cerr << "  [WALL] invert_wall=" << _wallSec << "s (NTL_time=" << (t1-t0) << "s)" << std::endl;
    std::cerr << "  [A-DEC TIME] rnd=" << rnd << " invert=" << t1-t0
               << "s bitCap=" << eData[0].bitCapacity() << std::endl;
    std::cerr << "  [noise] after invert(dec): " << eData[0].getNoiseBound()
               << " log2=" << log(eData[0].getNoiseBound())/log(2.0) << std::endl;

    for (long k=0; k<(long)eData.size(); k++) eData[k] += aesKey[rnd][k];
    std::cerr << "  [A-DEC TIME] rnd=" << rnd << " after addkey bitCap="
               << eData[0].bitCapacity() << std::endl;
  }
  std::cerr << "[A-DEC TOTAL] unpack=" << t_unpack << "s combinedLinearDecrypt=" << t_linear
             << "s pack=" << t_pack << "s invert=" << t_invert
             << "s boot=" << t_boot << "s bootCount=" << bootCount << std::endl;
  std::cerr << "[A-DEC TOTAL][WALL] invert_wall=" << wall_invert
             << "s unpack_wall=" << wall_unpack << "s linear_wall=" << wall_linear
             << "s pack_wall=" << wall_pack << "s boot_wall=" << wall_boot << "s" << std::endl;
}
