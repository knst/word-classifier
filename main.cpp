#include <fstream>
#include <iostream>

#include <stdexcept>

#include <cmath>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

const size_t bloomSize = 459983;
const size_t bloomK = 1;
const size_t bloomPrefixSize = 7;

vector<float> bigramProbability;

size_t bigramIndex(const string& bigram) {
    return (bigram[0] - 'a') * 26 + bigram[1] - 'a';
}

uint64_t hashString(const string& word, size_t index) {
    const uint64_t constants[] = {27644437, 115249, 33391, 108301, 115249};
    const size_t constantsCount = sizeof(constants) / sizeof(constants[0]);
    if (index + 1 > constantsCount)
        throw std::invalid_argument("overflow in hashString: index too long");

    uint64_t result = 0;
    for (auto i : word) {
        result = (result + i) * constants[index] + constants[index + 1];
    }
    return result % bloomSize;
}

void setBloom(const string& word, vector<bool>& bloom) {
    string prefix = word;
    if (prefix.size() > bloomPrefixSize)
        prefix = prefix.substr(0, bloomPrefixSize);
    for (size_t k = 0; k < bloomK; ++k) {
        bloom[hashString(prefix, k)] = true;
    }
}
bool testBloom(const string& word, const vector<bool>& bloom) {
    string prefix = word;
    if (prefix.size() > bloomPrefixSize)
        prefix = prefix.substr(0, bloomPrefixSize);
    for (size_t k = 0; k < bloomK; ++k) {
        if (!bloom[hashString(prefix, k)])
            return false;
    }
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

size_t countPair(const string& word) {
    size_t n = word.size();
    size_t count = 0;
    for (size_t i = 1; i < n; ++i) {
        if (word[i] == word[i - 1])
            ++count;
    }
    return count;
}

bool testWord(const string& word, const vector<bool>& bloom, bool isWord) {
    size_t n = word.size();
    string wordForBloom = word;
    bool isSword = false;

    if (n > 2 && word[n - 1] == 's' && word[n-2] == '\'') {
        wordForBloom = word.substr(0, n - 2);
        isSword = true;
    }
    if (n < 3)
        return true;
    if (n < 4 && isSword)
        return true;
    if (n > 2 && !testBloom(wordForBloom, bloom))
        return false;
//    if (isSword && !testWord(wordForBloom, bloom, isWord))
 //       return false;

    size_t m = isSword
        ? n - 2
        : n;
    size_t containAp = countAp(word);
    size_t strangeCount = 0;
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
    if (strangeCount > 1)
        return false;

//    if (countPair(word) > 2)
//        return false;

    float bigramProbSum = 0.0;
    float bigramProb = 1.0;
    float bigramSqrt = 0.0;
    for (size_t i = 1; i < m; ++i) {
        float probability = bigramProbability[bigramIndex(word.substr(i - 1, 2))];
        if (probability < 0.0000030)
            return false;
        bigramProbSum += probability;
        bigramProb *= probability;
        bigramSqrt += sqrt(probability);
    }
    bigramProb = log(bigramProb);
    float bigramProbMax[] = {
        0,
        0,
        0,
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
        -87, // 15
    };
    if (m >= 3 && bigramProb < bigramProbMax[m]) {
        return false;
    }
    float bigramSumMin[] = {
        0.008, // 10
        0.012, // 10
        0.015,
        0.024,
        0.035,
        0.043,
        0.045,
        0.050,
        0.060,
        0.070,
        0.080,
        0.090, // 20
        0.100,
    };
    if (m >= 9) {
        if (bigramProbSum *.69< bigramSumMin[m - 9]) {
            return false;
        }
    }
    if (m >= 15) {
        if (bigramProbSum > 0.17) {
            return false;
        }
    }

    float bigramSumSqrt[] = {
        0,
        0,
        0,
        0,
        0,
        0, // 5
        0.04,
        0.07,
        0.1,
        0.2,
        0.3, // 10
        0.4,
        0.6,
        0.8,
        0.9,
        1.1, // 15
        1.2,
        1.6,
        1.8,
        2.0, // 19
        2.2,
    };

    if (bigramSqrt *1.05 < bigramSumSqrt[m]) {
//    if (bigramSqrt *.75 < bigramSumSqrt[m]) {
        return false;
    }
/*
    vector<float> changedBigramProb;
    for (size_t i = 0; i < 2; ++i) {
        string wordC = word;
        for (size_t c = 'a'; c <= 'z'; ++c) {
            wordC[i] = c;

            float bigramProb = 0.0;
            for (size_t i = 1; i < m; ++i) {
                float probability = bigramProbability[bigramIndex(wordC.substr(i - 1, 2))];
                bigramProb += probability;
            }
//            bigramProb = log(bigramProb);
            changedBigramProb.push_back(bigramProb);
        }
    }
    sort(changedBigramProb.begin(), changedBigramProb.end());
    size_t index = 0;
    for (size_t i = 0; i < changedBigramProb.size(); ++i) {
        if (changedBigramProb[i] < bigramProbSum)
            ++index;
    }
//    if (static_cast<float>(index) / changedBigramProb.size() > 0.96)
//    if (static_cast<float>(index) / changedBigramProb.size() < 0.30)
//        return false;
*/
    return true;
}

void addBigram(const string& word, map<string, size_t>& biMap) {
    size_t n = word.size();
    for (size_t i = 1; i < n; ++i) {
        string bi;
        bi += word[i - 1];
        bi += word[i];
        ++biMap[bi];
    }
}

vector<float> countPropabilityBigram(map<string, size_t> bigramCount) {
    size_t summary = 0;
    for (auto& i : bigramCount) {
        if (i.second > 65535) {
            cerr << i.second << endl;
            i.second = 65535;
        }
        summary += i.second;
    }
    summary += 40000;
    ofstream probabilitiesFile("bigrams.bin", ios::binary);
    for (char c1 = 'a'; c1 <= 'z'; ++c1) {
        for (char c2 = 'a'; c2 <= 'z'; ++c2) {
            string bi;
            bi += c1;
            bi += c2;
            uint16_t count = static_cast<uint16_t>(bigramCount[bi]);
            probabilitiesFile.write(reinterpret_cast<char*>(&count), 2);
        }
    }

    vector<float> result(26 * 26);
    for (auto i : bigramCount) {
        result[bigramIndex(i.first)] = 1.0 * i.second / summary;
    }
    return result;
}

void printBiStat(const string& word, bool isWord) {
    float bigramProbSum = 0.0;
    float bigramProb = 1.0;
    float bigramSquare = 0.0;
    float bigramSqrt = 0.0;
    size_t n = word.size();
    if (n < 3)
        return ;
    if (n > 25)
        return ;
    for (size_t i = 1; i < n; ++i) {
        float probability = bigramProbability[bigramIndex(word.substr(i - 1, 2))];
        bigramProbSum += probability;
        bigramProb *= probability;
        bigramSquare += probability * probability;
        bigramSqrt += sqrt(probability);
    }
    bigramProb = log(bigramProb);
//    bigramSquare = log(bigramSquare);
//    bigramSqrt = log(bigramSqrt);
    cout << "sum: " << isWord << " " << bigramProbSum << ' ' << n << '\n';
    cout << "prob: " << isWord << " " << bigramProb << ' ' << n << '\n';
    cout << "square: " << isWord << " " << bigramSquare << ' ' << n << '\n';
    cout << "sqrt: " << isWord << " " << bigramSqrt << ' ' << n << '\n';
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "usage: " << argv[0] << " DIR" << endl;
        return 1;
    }

    size_t countTrashOk = 0;
    size_t countTrashFail = 0;
    size_t countWordsOk = 0;
    size_t countWordsFail = 0;

    string lex;

    string dictionaryFilename = "words-3.txt";
    ifstream dictionaryFile(dictionaryFilename.c_str());
    set<string> dictionary;

    map<string, size_t> wordsBi;

    while (dictionaryFile >> lex) {
        dictionary.insert(lex);
        addBigram(lex, wordsBi);
    }

    bigramProbability = countPropabilityBigram(wordsBi);

    vector<bool> filledBloom(bloomSize, true);
    vector<bool> bloom(bloomSize);
    size_t ok = 0;
    size_t fail = 0;
    for (auto i : dictionary) {
        if (testWord(i, filledBloom, true)) {
            ++ok;
            setBloom(i, bloom);
        } else {
            ++fail;
        }
    }
    cerr << ok << ' ' << fail << endl;
    vector<uint8_t> values;
    for (size_t i = 0; i < bloomSize; i += 8) {
        uint8_t value = 0;
        for (size_t j = 0; j < 8; ++j) {
            value = value << 1;
            if (i + j < bloomSize)
                value = value | bloom[i + j];
        }
        values.push_back(value);
    }
    ofstream bloomFilter("bloom.bin");
    for (auto i : values)
        bloomFilter << i;
    string dir = argv[1];
    string trashFile = dir + "/trash.list";
    string wordsFile = dir + "/words.list";
    ifstream trash(trashFile.c_str());
    ifstream words(wordsFile.c_str());

    while (trash >> lex) {
//        printBiStat(lex, false);
        if (testWord(lex, bloom, false)) {
            ++countTrashFail;
            if (rand() < 200000)
                cerr << lex << endl;
        } else {
            ++countTrashOk;
        }
    }
    cerr << "----" << endl;
    while (words >> lex) {
//        printBiStat(lex, true);
        if (testWord(lex, bloom, true)) {
            ++countWordsOk;
        } else {
            ++countWordsFail;
            if (rand() < 800000)
                cerr << lex << endl;
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
