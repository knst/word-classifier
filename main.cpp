#include <fstream>
#include <iostream>

#include <stdexcept>

#include <set>
#include <string>
#include <vector>

using namespace std;

const size_t bloomSize = 480000;
const size_t bloomK = 1;

uint64_t hashString(const string& word, size_t index) {
    const uint64_t constants[] = {27644437, 115249, 33391, 108301, 115249};
    const size_t constantsCount = sizeof(constants) / sizeof(constants[0]);
    if (index + 1 > constantsCount)
        throw std::invalid_argument("overflow in hashString: index too long");

    uint64_t result = 0;
    for (auto i : word) {
        result = result * constants[index] + constants[index + 1] + i;
    }
    return result % bloomSize;
}

void setBloom(const string& word, vector<bool>& bloom) {
    for (size_t k = 0; k < bloomK; ++k) {
        bloom[hashString(word, k)] = true;
    }
}
bool testBloom(const string& word, const vector<bool>& bloom) {
    for (size_t k = 0; k < bloomK; ++k) {
        if (!bloom[hashString(word, k)])
            return false;
    }
    return true;
}

bool testVowel(char i) {
    if (i == 'a' || i == 'e' || i == 'i' || i == 'o' || i == 'u')
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

bool testAp(const string& word) {
    for (auto i : word)
        if (i == '\'')
            return true;
    return false;
}

bool testTripple(const string& word) {
    size_t n = word.size();
    for (size_t i = 2; i < n; ++i) {
        if (word[i] == word[i - 1] && word[i] == word[i - 2])
            return true;
    }
    return false;
}

bool testTrippleVowel(const string& word) {
    size_t n = word.size();
    for (size_t i = 2; i < n; ++i) {
        if (testVowel(word[i]) && testVowel(word[i-1]) && testVowel(word[i-2]))
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

size_t endedByEd(const string& word) {
    size_t n = word.size();
    if (n < 4)
        return false;
    if (word[n-3] == 'd' && word[n-4] == 'e')
        return true;
    return false;
}

bool testWord(const string& word, const vector<bool>& bloom) {
    if (!testBloom(word, bloom))
        return false;

    size_t n = word.size();
    bool containAp = testAp(word);
    if (n < 3 && containAp)
        return false;
    if (n < 3)
        return true;
    if (n == 3)
        return false;

    if (testTripple(word))
        return false;
//    if (testTrippleVowel(word) && n < 8)
//        return false;
    if (n > 14)
        return false;

    if (word[n - 1] == 's' && word[n - 2] == '\'') {
        return true;
    }

    size_t consolantCount = countConsolant(word);
    size_t vowelCount = n - consolantCount;
//    if (vowelCount == 1 && endedByEd(word))
//        return false;
    if (vowelCount == 0)
        return false;
    if (consolantCount == 0)
        return false;
    if (n > 7 && vowelCount < n / 5)
        return false;
    if (n > 7 && consolantCount < n / 4)
        return false;

//    if (word[0] >= 'z')
//        return false;

//    if (countPair(word) > 3)
 //       return false;
    if (containAp)
        return false;
    return true;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "usage: " << argv[0] << "DIR" << endl;
        return 1;
    }

    size_t countTrashOk = 0;
    size_t countTrashFail = 0;
    size_t countWordsOk = 0;
    size_t countWordsFail = 0;

    string lex;

    string dictionaryFilename = "words-uniq.txt";
    ifstream dictionaryFile(dictionaryFilename.c_str());
    set<string> dictionary;
    while (dictionaryFile >> lex) {
        dictionary.insert(lex);
    }

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
        if (testWord(lex, bloom)) {
            ++countTrashFail;
            if (rand() < 200000)
                cerr << lex << endl;
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

    double precisionTrash = 1.0 * countTrashOk / (countTrashOk + countTrashFail);
    double precisionWords = 1.0 * countWordsOk / (countWordsOk + countWordsFail);
    double precision = 1.0 * (countTrashOk + countWordsOk) / (countTrashOk + countTrashFail + countWordsFail + countWordsOk);
    cout << "Trash: " << precisionTrash << " " << countTrashOk << " : " << countTrashFail << endl;
    cout << "Words: " << precisionWords << " " << countWordsOk << " : " << countWordsFail << endl;
    cout << "Precision: " << precision << endl;
}
