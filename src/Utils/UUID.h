#ifndef UUID_H
#define UUID_H

#include <cstdint>

namespace AnasenSim {

    //Use szudzik pairing method to make unqiue key out of two unsigned ints.
	static constexpr uint32_t GetUUID(uint32_t a, uint32_t b)
	{
		return a >= b ? a*a + a + b : b*b + a;
	}

}

#endif