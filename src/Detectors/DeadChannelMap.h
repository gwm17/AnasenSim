#ifndef DEAD_CHANNEL_MAP_H
#define DEAD_CHANNEL_MAP_H

#include <cstdint>
#include <unordered_map>
#include <string>

namespace AnasenSim {

    class DeadChannelMap
    {
    public:

        enum class DetectorType
        {
            Barrel1,
            Barrel2,
            PC,
            None
        };

        enum class ChannelType
        {
            Front,
            Back,
            Wire,
            None
        };

        DeadChannelMap();
        ~DeadChannelMap();

        void ReadFile(const std::string& filename);

        //Return true if wire is dead, false if wire is alive
        bool IsWireDead(uint32_t wireid) const;
        //Methods for Si
        //Return true if channel is dead, false if channel is alive
        bool IsChannelDead(uint32_t detid, uint32_t channelid, DetectorType dettype, ChannelType chantype) const;
        //Return true if channel pair is dead, false if channel pair is alive
        bool IsChannelPairDead(uint32_t detid, uint32_t channelfront, uint32_t channelback, DetectorType dettype) const;

        bool IsValid() const { return m_isValid; }

    private:
        std::unordered_map<uint32_t, bool> m_barrel1Map;
        std::unordered_map<uint32_t, bool> m_barrel2Map;
        std::unordered_map<uint32_t, bool> m_pcMap;
        bool m_isValid;

        static constexpr uint32_t s_sx3BackChannelOffset = 4;
        static constexpr uint32_t s_pcDetID = 99;
    };

}

#endif