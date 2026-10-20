#pragma once
#include <stdint.h>
#include <helib/helib.h>

#ifdef USE_ZZX_POLY
#define PolyType NTL::ZZX
#else
#define PolyType helib::DoubleCRT
#endif

class HomAES {
  const helib::EncryptedArrayDerived<helib::PA_GF2> ea2;
  std::vector<PolyType> encAffMat, decAffMat;
  PolyType affVec;
  std::vector<PolyType> encLinTran, decLinTran;
  NTL::GF2X XinSlots;
  NTL::Mat<NTL::GF2X> unpacking;
public:
  static const NTL::GF2X aesPoly;
  explicit HomAES(const helib::Context& context);
  void batchRecrypt(std::vector<helib::Ctxt>& data) const;
  void setPackingConstants();
  void debugUnpackTest(const helib::PubKey& pk, const helib::SecKey& sk) const;  // 临时验证方法
  void batchRecryptIsolationTest(const helib::PubKey& pk, const helib::SecKey& sk) const;  // 隔离测试: 只测batchRecrypt是否保值
  void bitLevelRoundTestA(std::vector<helib::Ctxt>& eData,
                           const std::vector<std::vector<helib::Ctxt>>& aesKey) const;  // 方案A完整10轮测试
  void bitLevelRoundTestADecrypt(std::vector<helib::Ctxt>& eData,
                                  const std::vector<std::vector<helib::Ctxt>>& aesKey) const;  // 方案A完整10轮解密测试
  void testSboxOnly(std::vector<helib::Ctxt>& eData) const;  // 只测invert+affine是否等于标准S-box
  void encryptAESkey(std::vector<helib::Ctxt>& eKey, NTL::Vec<uint8_t>& aesKey, const helib::PubKey& hePK) const;
  void encryptAESkey_singleBlock16Ctxt(std::vector<std::vector<helib::Ctxt>>& eKeySchedule, NTL::Vec<uint8_t>& aesKey, const helib::PubKey& hePK) const;  // 方案A专用: 每轮密钥16字节各占独立密文slot0
  void encryptAESkey_broadcastBlock16Ctxt(std::vector<std::vector<helib::Ctxt>>& eKeySchedule, NTL::Vec<uint8_t>& aesKey, const helib::PubKey& hePK) const;  // 方案B专用: 每轮密钥16字节各广播复制到全部slot
  void encryptAESkey_multiGroupBroadcastCtxt(std::vector<std::vector<helib::Ctxt>>& eKeySchedule, NTL::Vec<uint8_t>& aesKey, long numGroups, const helib::PubKey& hePK) const;  // 方案C专用: numGroups组, 每组16密文, 密钥在组内广播且组间重复
  void encryptAESkey_bitLevel(std::vector<std::vector<std::vector<helib::Ctxt>>>& eKeyBits, NTL::Vec<uint8_t>& aesKey, const helib::PubKey& hePK) const;  // bit态轮密钥: [轮次][字节位置][bit位置], 每密文所有slot填同一值
  void homAESenc(std::vector<helib::Ctxt>& eData, const std::vector<helib::Ctxt>& eKey, const NTL::Vec<uint8_t> inBytes) const;
  void homAESdec(std::vector<helib::Ctxt>& eData, const std::vector<helib::Ctxt>& eKey, const NTL::Vec<uint8_t> inBytes) const;
  void homAESenc(std::vector<helib::Ctxt>& eData, const std::vector<helib::Ctxt>& eKey) const;
  void homAESdec(std::vector<helib::Ctxt>& eData, const std::vector<helib::Ctxt>& eKey) const;
  const helib::EncryptedArrayDerived<helib::PA_GF2>& getEA() const { return ea2; }
};

void encode4AES(NTL::Vec<NTL::ZZX>& encData, const NTL::Vec<uint8_t>& data, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES(NTL::Vec<uint8_t>& data, const NTL::Vec<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_singleBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t block[16], const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES_singleBlock16Ctxt(uint8_t block[16], const std::vector<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_singleBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t block[16], const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES_singleBlock16Ctxt(uint8_t block[16], const std::vector<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_singleBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t block[16], const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES_singleBlock16Ctxt(uint8_t block[16], const std::vector<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_multiBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t* blocks, long numBlocks, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES_multiBlock16Ctxt(uint8_t* blocks, long numBlocks, const std::vector<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_multiGroupBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t* blocks, long numGroups, long numBlocksPerGroup, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void decode4AES_multiGroupBlock16Ctxt(uint8_t* blocks, long numGroups, long numBlocksPerGroup, const std::vector<NTL::ZZX>& encData, const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);
void encode4AES_broadcastBlock16Ctxt(std::vector<NTL::ZZX>& encData, const uint8_t block[16], const helib::EncryptedArrayDerived<helib::PA_GF2>& ea2);


