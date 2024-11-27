#include "coin.h"
#include "util.h"
#include "amount.h"

#include <openssl/rand.h>
#include <sstream>
#include "openssl_context.h"
#include "primitives/mint_spend.h"

namespace sigma {

std::ostream& operator<<(std::ostream& stream, CoinDenomination denomination) {
    int64_t denom_value;
    DenominationToInteger(denomination, denom_value);
    stream << denom_value;
    return stream;
}

bool DenominationToInteger(CoinDenomination denom, int64_t& denom_out) {
    CValidationState dummy_state;
    return DenominationToInteger(denom, denom_out, dummy_state);
}

bool DenominationToInteger(CoinDenomination denom, int64_t& denom_out, CValidationState &state) {

    return false;
}

bool RealNumberToDenomination(const double& value, CoinDenomination& denom_out) {
    return IntegerToDenomination(value * COIN, denom_out);
}

bool StringToDenomination(const std::string& str, CoinDenomination& denom_out) {

    return false;
}

std::string DenominationToString(const CoinDenomination& denom) {

    throw std::runtime_error("Unsupported denomination, unable to convert to string.");
}

bool IntegerToDenomination(int64_t value, CoinDenomination& denom_out) {
    CValidationState dummy_state;
    return IntegerToDenomination(value, denom_out, dummy_state);
}

bool IntegerToDenomination(int64_t value, CoinDenomination& denom_out, CValidationState &state) {

    return false;
}

void GetAllDenoms(std::vector<sigma::CoinDenomination>& denominations_out) {
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_100);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_25);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_10);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_1);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_0_5);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_0_1);
    denominations_out.push_back(CoinDenomination::SIGMA_DENOM_0_05);
}

//class PublicCoin
PublicCoin::PublicCoin()
    : denomination(CoinDenomination::SIGMA_DENOM_1)
{

}

PublicCoin::PublicCoin(const GroupElement& coin, const CoinDenomination d)
    : value(coin)
    , denomination(d)
{
}

const GroupElement& PublicCoin::getValue() const{
    return this->value;
}

uint256 const & PublicCoin::getValueHash() const {
    if(valueHash.IsNull())
        valueHash = primitives::GetPubCoinValueHash(value);
    return valueHash;
}


CoinDenomination PublicCoin::getDenomination() const {
    return denomination;
}

bool PublicCoin::operator==(const PublicCoin& other) const{
    return (*this).value == other.value;
}

bool PublicCoin::operator!=(const PublicCoin& other) const{
    return (*this).value != other.value;
}

bool PublicCoin::validate() const{
    return this->value.isMember() && !this->value.isInfinity();
}

//class PrivateCoin
PrivateCoin::PrivateCoin(const Params* p, CoinDenomination denomination, int version)
    : params(p)
{
        this->version = version;
        this->mintCoin(denomination);
}

const Params * PrivateCoin::getParams() const {
    return this->params;
}

const PublicCoin& PrivateCoin::getPublicCoin() const{
    return this->publicCoin;
}

const Scalar& PrivateCoin::getSerialNumber() const{
    return this->serialNumber;
}

const Scalar& PrivateCoin::getRandomness() const{
    return this->randomness;
}

const unsigned char* PrivateCoin::getEcdsaSeckey() const {
     return this->ecdsaSeckey;
}

void PrivateCoin::setEcdsaSeckey(const std::vector<unsigned char> &seckey) {
    if (seckey.size() == sizeof(ecdsaSeckey))
        std::copy(seckey.cbegin(), seckey.cend(), &ecdsaSeckey[0]);
    else
        throw std::invalid_argument("EcdsaSeckey size does not match.");
}

void PrivateCoin::setEcdsaSeckey(uint256 &seckey) {
    if (seckey.size() == sizeof(ecdsaSeckey))
        std::copy(seckey.begin(), seckey.end(), &ecdsaSeckey[0]);
    else
        throw std::invalid_argument("EcdsaSeckey size does not match.");
}

unsigned int PrivateCoin::getVersion() const {
    return this->version;
}

void PrivateCoin::setPublicCoin(const PublicCoin& p) {
    publicCoin = p;
}

void PrivateCoin::setRandomness(const Scalar& n) {
    randomness = n;
}

void PrivateCoin::setSerialNumber(const Scalar& n) {
    serialNumber = n;
}

void PrivateCoin::setVersion(unsigned int nVersion){
    version = nVersion;
}

void PrivateCoin::mintCoin(const CoinDenomination denomination){
    // Create a key pair
    secp256k1_pubkey pubkey;
    do {
        if (RAND_bytes(this->ecdsaSeckey, sizeof(this->ecdsaSeckey)) != 1) {
            throw std::runtime_error("Unable to generate randomness");
        }
    } while (!secp256k1_ec_pubkey_create(
        OpenSSLContext::get_context(), &pubkey, this->ecdsaSeckey));

    // Hash the public key in the group to obtain a serial number
    serialNumber = serialNumberFromSerializedPublicKey(
        OpenSSLContext::get_context(), &pubkey);

    randomness.randomize();
    GroupElement commit = SigmaPrimitives<Scalar, GroupElement>::commit(
            params->get_g(), serialNumber, params->get_h0(), randomness);
    publicCoin = PublicCoin(commit, denomination);
}

Scalar PrivateCoin::serialNumberFromSerializedPublicKey(
        const secp256k1_context *context,
        secp256k1_pubkey *pubkey) {
    std::vector<unsigned char> pubkey_hash(32, 0);

    static const unsigned char one[32] = {
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01
    };

    // We use secp256k1_ecdh instead of secp256k1_serialize_pubkey to avoid a timing channel.
    if (1 != secp256k1_ecdh(context, pubkey_hash.data(), pubkey, &one[0])) {
        throw std::runtime_error("Unable to compute public key hash with secp256k1_ecdh.");
    }

	std::string zpts(PRIVCOIN_PUBLICKEY_TO_SERIALNUMBER);
	std::vector<unsigned char> pre(zpts.begin(), zpts.end());
    std::copy(pubkey_hash.begin(), pubkey_hash.end(), std::back_inserter(pre));

	unsigned char hash[CSHA256::OUTPUT_SIZE];
    CSHA256().Write(pre.data(), pre.size()).Finalize(hash);

    // Use 32 bytes of hash as coin serial.
    return Scalar(hash);
}

} // namespace sigma

namespace std {

string to_string(::sigma::CoinDenomination denom)
{
    throw invalid_argument("the specified denomination is not valid");
}

} // namespace std
