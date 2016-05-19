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

const size_t bloomSize = 460000;
const size_t bloomK = 1;

map<char, double> letterProbability;
vector<double> bigramProbability;

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
    if (prefix.size() > 7)
        prefix = prefix.substr(0, 7);
    for (size_t k = 0; k < bloomK; ++k) {
        bloom[hashString(prefix, k)] = true;
    }
}
bool testBloom(const string& word, const vector<bool>& bloom) {
    string prefix = word;
    if (prefix.size() > 7)
        prefix = prefix.substr(0, 7);
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

size_t countConsolant(const string& word) {
    size_t count = 0;
    for (auto i : word) {
        if (testVowel(i))
            ++count;
    }
    return count;
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
    if (isSword && !testWord(wordForBloom, bloom, isWord))
        return false;

    size_t m = isSword
        ? n - 2
        : n;
    if (m <= 6)
        return true;
    size_t containAp = countAp(word);
    size_t strangeCount = 0;
    if (containAp && !isSword)
        return false;
    if (containAp > 1)
        return false;

    if (testTripple(word) && n != 3)
        return false;


    if (n > 25) {
        // too many correct words is lost, but trash filtered better (at this moment).
        return false;
    }

    size_t consolantSeq = countConsolantSeq(word);
    size_t vowelSeq = countVowelSeq(word);
    size_t consolantCount = countConsolant(word);
    size_t vowelCount = n - consolantCount;
    if (consolantSeq > 4)
        return false;
    if (vowelSeq > 6)
        return false;
    if (n > 7 && vowelCount < n / 6)
        ++strangeCount;
    if (n > 7 && consolantCount < n / 5)
        ++strangeCount;
    if (strangeCount > 1)
        return false;

   // if (word[0] >= 'z')
    //    ++strangeCount;

    if (countPair(word) > 2)
        return false;

//    if (containAp)
//        return false;
    double bigramProbSum = 0.0;
    double bigramProb = 1.0;
    double bigramSqrt = 0.0;
    for (size_t i = 1; i < m; ++i) {
        double probability = bigramProbability[bigramIndex(word.substr(i - 1, 2))];
        if (probability < 0.0000030)
            return false;
        bigramProbSum += probability;
        bigramProb *= probability;
        bigramSqrt += sqrt(probability);
    }
    bigramProb = log(bigramProb);
    double bigramProbMax[] = {
        0,
        0,
        0,
        0,
        -20,
        -30, // 5
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
    if (m >= 7 && bigramProb < bigramProbMax[m]) {
        return false;
    }
    double bigramSumMin[] = {
        0.006, // 10
        0.012, // 10
        0.015,
        0.020,
        0.030,
        0.035,
        0.045,
        0.050,
        0.060,
        0.070,
        0.080,
        0.090, // 20
        0.100,
    };
    if (m >= 9) {
        if (bigramProbSum *.75< bigramSumMin[m - 9]) {
            return false;
        }
    }
    if (m >= 15) {
        if (bigramProbSum > 0.17) {
            return false;
        }
    }

    double bigramSumSqrt[] = {
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

    double letterProbSum = 0.0;
    double letterProb = 1.0;
    for (auto i : word) {
        double probability = letterProbability[i];
        letterProbSum += probability;
        letterProb *= probability;
    }
    letterProb = log(letterProb);
    double letterSum[] = {
        0,
        0,
        0,
        0, // 3
        0.05,
        0.1, // 5
        0.15,
        0.20,
        0.25,
        0.30,
        0.35, // 10
        0.40,
        0.45,
        0.50,
        0.55,
        0.60,
        0.65,
        0.70,
        0.75,
        0.80,
        0.85, // 20
    };
    if (letterProbSum < letterSum[m]) {
//        return false;
    }
/*
    vector<double> changedBigramProb;
    for (size_t i = 0; i < m; ++i) {
        string wordC = word;
        for (size_t c = 'a'; c <= 'z'; ++c) {
            wordC[i] = c;

            double bigramProb = 0.0;
            for (size_t i = 1; i < m; ++i) {
                double probability = bigramProbability[bigramIndex(wordC.substr(i - 1, 2))];
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
    cout << "proprtion: " << isWord << ' ' << static_cast<double>(index) / changedBigramProb.size() << ' ' << n << '\n';
    */
//    if (changedSumBigramProb / changedCount < bigramProb * 0.9) {
//        return false;
//    }

    return true;
}

void initBigram(map<string, size_t>& biMap) {
    for (char c = 'a'; c <= 'z'; ++c) {
        for (char c2 = 'a'; c2 <= 'z'; ++c2) {
            string bi;
            bi += c;
            bi += c2;
            biMap[bi] = 0;
        }
    }
}

void addLetter(const string& word, map<char, size_t>& letterMap) {
    for (auto i : word)
        ++letterMap[i];
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

map<char, double> countPropabilityLetter(map<char, size_t> letterCount) {
    size_t summary = 0;
    for (auto i : letterCount)
        summary += i.second;
    map<char, double> result;
    for (auto i : letterCount) {
        result[i.first] = 1.0 * i.second / summary;
    }
    return result;
}

vector<double> countPropabilityBigram(map<string, size_t> bigramCount) {
    size_t summary = 0;
    for (auto i : bigramCount)
        summary += i.second;
    vector<double> result(26 * 26);
    for (auto i : bigramCount) {
        result[bigramIndex(i.first)] = 1.0 * i.second / summary;
    }
    return result;
}

void printBiStat(const string& word, bool isWord) {
    double bigramProbSum = 0.0;
    double bigramProb = 1.0;
    double bigramSquare = 0.0;
    double bigramSqrt = 0.0;
    size_t n = word.size();
    if (n < 3)
        return ;
    if (n > 25)
        return ;
    for (size_t i = 1; i < n; ++i) {
        double probability = bigramProbability[bigramIndex(word.substr(i - 1, 2))];
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

void printLetterStat(const string& word, bool isWord) {
    size_t n = word.size();
    if (n > 25)
        return ;
    double letterProbSum = 0.0;
    double letterProb = 1.0;
    for (auto i : word) {
        double probability = letterProbability[i];
        letterProbSum += probability;
        letterProb *= probability;
    }
    letterProb = log(letterProb);
    cout << "sum-l: " << isWord << " " << letterProbSum << ' ' << n << '\n';
    cout << "prob-l: " << isWord << " " << letterProb << ' ' << n << '\n';
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
    map<string, size_t> trashBi;
    map<char, size_t> letter;
    initBigram(wordsBi);
    initBigram(trashBi);

    while (dictionaryFile >> lex) {
        addLetter(lex, letter);
        addBigram(lex, wordsBi);
        dictionary.insert(lex);
    }

    letterProbability = countPropabilityLetter(letter);
    bigramProbability = countPropabilityBigram(wordsBi);

    vector<bool> bloom(bloomSize);
    for (auto i : dictionary) {
        setBloom(i, bloom);
    }

    string dir = argv[1];
    string trashFile = dir + "/trash.list";
    string wordsFile = dir + "/words.list";
    ifstream trash(trashFile.c_str());
    ifstream words(wordsFile.c_str());

    while (trash >> lex) {
//        printBiStat(lex, false);
//        printLetterStat(lex, false);
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
 //       printLetterStat(lex, true);
        if (testWord(lex, bloom, true)) {
            ++countWordsOk;
        } else {
            ++countWordsFail;
            if (rand() < 800000)
                cerr << lex << endl;
        }
    }

    double precisionTrash = 1.0 * countTrashOk / (countTrashOk + countTrashFail);
    double precisionWords = 1.0 * countWordsOk / (countWordsOk + countWordsFail);
    double precision = 1.0 * (countTrashOk + countWordsOk) / (countTrashOk + countTrashFail + countWordsFail + countWordsOk);
    cout << "Trash: " << precisionTrash << " " << countTrashOk << " : " << countTrashFail << endl;
    cout << "Words: " << precisionWords << " " << countWordsOk << " : " << countWordsFail << endl;
    cout << "Precision: " << precision << endl;
    return 0;
}
