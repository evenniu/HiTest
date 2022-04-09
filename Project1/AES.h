#pragma once
#include <Eigen/Dense>

#define DLLEXPORT __declspec(dllexport)

namespace MyApp
{
	class DLLEXPORT  AES
	{
	private:
		unsigned char m_Sbox[256];
		unsigned char m_InvSbox[256];
		unsigned char m_w[11][4][4];
	public:
		AES(void);
		AES(unsigned char* key);
		unsigned char* Cipher(unsigned char* input);
		void AddRoundKey(unsigned char state[][4], unsigned char k[][4]);
		~AES();
		//unsigned char* Cipher(unsigned char* input);
		void KeyExpansion(unsigned char* key, unsigned char w[][4][4]);
		int strToUChar(const char* ch, unsigned char* uch);
		int add(int a, int b);
		int TestEigen();
	};
}
