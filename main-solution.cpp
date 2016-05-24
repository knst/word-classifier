#include <fstream>
#include <iostream>

#include <cmath>

#include <string>
#include <vector>

using namespace std;

const size_t bloomSize = 459983;
const size_t bloomPrefixSize = 7;
const size_t bigramCount = 26 * 26;

vector<float> bigramProbability;

size_t bigramIndex(const string& bigram) {
    return (bigram[0] - 'a') * 26 + bigram[1] - 'a';
}

uint64_t hashString(const string& word) {
    const uint64_t constants[] = {27644437, 115249, 33391, 108301, 115249};

    uint64_t result = 0;
    for (auto i : word) {
        result = (result + i) * constants[0] + constants[1];
    }
    return result % bloomSize;
}

bool testBloom(const string& word, const vector<bool>& bloom) {
    string prefix = word;
    if (prefix.size() > bloomPrefixSize)
        prefix = prefix.substr(0, bloomPrefixSize);
    if (!bloom[hashString(prefix)])
        return false;
    return true;
}

bool testVowel(char i) {
    if (i == 'a' || i == 'e' || i == 'i' || i == 'o' || i == 'u' || i =='y')
        return true;
    return false;
}

size_t countConsolantSeq(const string& word) {
    size_t bestCount = 0;
    size_t n = word.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = i;
        while (j < n && testVowel(word[j]))
            ++j;
        bestCount = max(bestCount, j - i);
    }
    return bestCount;
}
size_t countVowelSeq(const string& word) {
    size_t bestCount = 0;
    size_t n = word.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = i;
        while (j < n && !testVowel(word[j]))
            ++j;
        bestCount = max(bestCount, j - i);
    }
    return bestCount;
}

size_t countAp(const string& word) {
    size_t count = 0;
    for (auto i : word)
        if (i == '\'')
            ++count;
    return count;
}

bool testTripple(const string& word) {
    size_t n = word.size();
    for (size_t i = 2; i < n; ++i) {
        if (word[i] == word[i - 1] && word[i] == word[i - 2])
            return true;
    }
    return false;
}

bool testWord(const string& word, const vector<bool>& bloom) {
    size_t n = word.size();
    string wordForBloom = word;
    bool isSword = false;

    if (n > 2 && word[n - 1] == 's' && word[n-2] == '\'') {
        wordForBloom = word.substr(0, n - 2);
        isSword = true;
    }
    size_t m = isSword
        ? n - 2
        : n;
    size_t containAp = countAp(word);

    if (n < 3 && containAp)
        return false;
    if (n == 2 && bigramProbability[bigramIndex(word)] < 3.0e-6)
        return false;
    if (n < 3)
        return true;
    if (n < 4 && isSword)
        return true;
    if (n > 2 && !testBloom(wordForBloom, bloom))
        return false;

    if (containAp && !isSword)
        return false;
    if (containAp > 1)
        return false;

    if (testTripple(word) && n > 3)
        return false;

    size_t consolantSeq = countConsolantSeq(word);
    size_t vowelSeq = countVowelSeq(word);
    if (consolantSeq > 5)
        return false;
    if (vowelSeq > 5)
        return false;

    float bigramProbSum = 0.0;
    float bigramProb = 1.0;
    float bigramSqrt = 0.0;
    for (size_t i = 1; i < m; ++i) {
        float probability = bigramProbability[bigramIndex(word.substr(i - 1, 2))];
        if (probability < 3.0e-6)
            return false;
        bigramProbSum += probability;
        bigramProb *= probability;
        bigramSqrt += sqrt(probability);
    }
    bigramProb = log(bigramProb);
    if (n > 15) {
        return false;
    }
    float bigramProbMax[] = {
        -21,
        -25,
        -32, // 5
        -37,
        -45,
        -51,
        -62,
        -69, // 10
        -74,
        -77,
        -81,
        -84,
        -87, // 15,
        -89,
        -91,
    };
    if (bigramProb < bigramProbMax[n - 3]) {
        return false;
    }
    float bigramSumMin[] = {
        0.005,
        0.008, // 10
        0.012, // 10
        0.015,
        0.024,
        0.035,
        0.043,
        0.045,
        0.050,
    };
    if (m >= 8) {
        if (bigramProbSum *.69 < bigramSumMin[m - 8]) {
            return false;
        }
    }

    float bigramSumSqrt[] = {
        0.29, // 9
        0.36, // 10
        0.42,
        0.57,
        0.7,
        0.86,
        1.05, // 15
    };

    if (m > 8 && bigramSqrt < bigramSumSqrt[m - 9]) {
        return false;
    }
    return true;
}

vector<float> countProbabilityBigram() {
    vector<float> result(bigramCount);
    ifstream probabilitiesFile("bigrams.bin", ios::binary);
    size_t summary = 40000;
    for (size_t i = 0; i < bigramCount; ++i) {
        uint16_t count;
        probabilitiesFile.read(reinterpret_cast<char *>(&count), 2);
        summary += count;
        result[i] = count;
    }

    for (size_t i = 0; i < bigramCount; ++i) {
        result[i] = 1. * result[i] / summary;
    }
    return result;
}

int main(int argc, char *argv[]) {
    bigramProbability = countProbabilityBigram();

    vector<bool> bloom(bloomSize + 7);
    vector<uint8_t> values;
    ifstream bloomFile("bloom.bin", ios::binary);
    for (size_t i = 0; i < bloomSize; i += 8) {
        uint8_t value = 0;
        bloomFile.read(reinterpret_cast<char *>(&value), 1);
        for (size_t j = 7; j < 8; --j) {
            bloom[i + j] = value & 1;
            value = value >> 1;
        }
    }

    if (argc < 2) {
        cerr << "usage: " << argv[0] << " DIR" << endl;
        return 1;
    }

    string dir = argv[1];
    string trashFile = dir + "/trash.list";
    string wordsFile = dir + "/words.list";
    ifstream trash(trashFile.c_str());
    ifstream words(wordsFile.c_str());

    size_t countTrashOk = 0;
    size_t countTrashFail = 0;
    size_t countWordsOk = 0;
    size_t countWordsFail = 0;

    string lex;
    while (trash >> lex) {
        if (testWord(lex, bloom)) {
            ++countTrashFail;
        } else {
            ++countTrashOk;
        }
    }
    while (words >> lex) {
        if (testWord(lex, bloom)) {
            ++countWordsOk;
        } else {
            ++countWordsFail;
        }
    }

    float precisionTrash = 100.0 * countTrashOk / (countTrashOk + countTrashFail);
    float precisionWords = 100.0 * countWordsOk / (countWordsOk + countWordsFail);
    float precision = 100.0 * (countTrashOk + countWordsOk) / (countTrashOk + countTrashFail + countWordsFail + countWordsOk);
    cout << "Trash: " << precisionTrash << " " << countTrashOk << " : " << countTrashFail << endl;
    cout << "Words: " << precisionWords << " " << countWordsOk << " : " << countWordsFail << endl;
    cout << "Precision: " << precision << endl;
    return 0;
}
