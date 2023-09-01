#include <iostream>

// Base class for different versions of Beta
class BaseBeta {
public:
    BaseBeta(int arg) : argument(arg) {}
    virtual void someFunction() = 0; // Pure virtual function
    // Other common members and functions for BaseBeta
protected:
    int argument;
};

// Different versions of Beta derived from BaseBeta
class BetaV1 : public BaseBeta {
public:
    BetaV1(int arg) : BaseBeta(arg) {}
    void someFunction() override {
        std::cout << "BetaV1::someFunction(), Argument: " << argument << std::endl;
    }
    // Other members and functions specific to BetaV1
};

class BetaV2 : public BaseBeta {
public:
    BetaV2(int arg) : BaseBeta(arg) {}
    void someFunction() override {
        std::cout << "BetaV2::someFunction(), Argument: " << argument << std::endl;
    }
    // Other members and functions specific to BetaV2
};

// Alpha class remains unchanged
class Alpha {
private:
    BaseBeta* betaInstance; // Pointer to BaseBeta, can hold any version of Beta

public:
    Alpha(BaseBeta* beta) : betaInstance(beta) {}

    void useBeta() {
        betaInstance->someFunction(); // Delegate the call to the wrapped Beta instance
    }
    // Other members and functions for Alpha
};

// Wrapper class WrapAlpha
class WrapAlpha {
private:
    Alpha alphaInstance;

public:
    WrapAlpha(const std::string& option, int arg) : alphaInstance(nullptr) {
        if (option == "V1") {
            alphaInstance = Alpha(new BetaV1(arg));
        } else if (option == "V2") {
            alphaInstance = Alpha(new BetaV2(arg));
        } else {
            std::cerr << "Invalid option: " << option << std::endl;
        }
    }

    void useAlpha() {
        alphaInstance.useBeta(); // Delegate the call to Alpha instance
    }
    // Other members and functions for WrapAlpha
};

int main() {
    WrapAlpha wrapAlphaV1("V1", 42); // Pass the argument here
    WrapAlpha wrapAlphaV2("V2", 73); // Pass the argument here

    wrapAlphaV1.useAlpha(); // Uses BetaV1's implementation
    wrapAlphaV2.useAlpha(); // Uses BetaV2's implementation

    return 0;
}
