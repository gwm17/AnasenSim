#include "DeadChannelMap.h"
#include "Utils/UUID.h"
#include <fstream>
#include <iostream>

namespace AnasenSim {

    DeadChannelMap::DetectorType ConvertStringToDetectorType(const std::string& type)
    {
        if(type == "barrel1")
            return DeadChannelMap::DetectorType::Barrel1;
        else if (type == "barrel2")
            return DeadChannelMap::DetectorType::Barrel2;
        else if (type == "pc")
            return DeadChannelMap::DetectorType::PC;
        else
            return DeadChannelMap::DetectorType::None;
    }

    DeadChannelMap::ChannelType ConvertStringToChannelType(const std::string& type)
    {
        if(type == "front")
            return DeadChannelMap::ChannelType::Front;
        else if (type == "back")
            return DeadChannelMap::ChannelType::Back;
        else if (type == "wire")
            return DeadChannelMap::ChannelType::Wire;
        else
            return DeadChannelMap::ChannelType::None;
    }

    DeadChannelMap::DeadChannelMap() :
        m_isValid(false)
    {
    }

    DeadChannelMap::~DeadChannelMap() {}

    void DeadChannelMap::ReadFile(const std::string& filename)
    {
        std::ifstream input(filename);
        if(!input.is_open())
        {
            m_isValid = false;
            return;
        }
        std::string junk;
        std::getline(input, junk);
        uint32_t detid, channel, uuid;
        DetectorType detType;
        ChannelType chanType;
        while(input >> detid)
        {
            input >> junk;
            detType = ConvertStringToDetectorType(junk);
            input >> junk;
            chanType = ConvertStringToChannelType(junk);
            input >> channel;
            if(detType == DetectorType::None)
            {
                std::cerr << "Error parsing dead channel map! Unidentified detector type" << std::endl;
                m_isValid = false;
                return;
            }
            else if(detType == DetectorType::Barrel1)
            {
                if(chanType == ChannelType::None)
                {
                    std::cerr << "Error parsing dead channel map! Unidentified channel type" << std::endl;
                    m_isValid = false;
                    return;
                }
                else if (chanType == ChannelType::Back)
                {
                    channel += s_sx3BackChannelOffset;
                }
                m_barrel1Map[GetUUID(detid, channel)] = true;
            }
            else if(detType == DetectorType::Barrel2)
            {
                if (chanType == ChannelType::None)
                {
                    std::cerr << "Error parsing dead channel map! Unidentified channel type" << std::endl;
                    m_isValid = false;
                    return;
                }
                else if (chanType == ChannelType::Back)
                {
                   channel += s_sx3BackChannelOffset;
                }
                m_barrel2Map[GetUUID(detid, channel)] = true;
            }
            else if(detType == DetectorType::PC)
            {
                if (chanType == ChannelType::None)
                {
                    std::cerr << "Error parsing dead channel map! Unidentified channel type" << std::endl;
                    m_isValid = false;
                    return;
                }
                m_pcMap[GetUUID(s_pcDetID, channel)] = true;
            }
        }

        input.close();
        m_isValid = true;
    }

    bool DeadChannelMap::IsWireDead(uint32_t wireID) const 
    {
        if(!m_isValid)
            return false;
        auto iter = m_pcMap.find(GetUUID(s_pcDetID, wireID));
        if (iter != m_pcMap.end())
            return true;
        return false;
    }

    bool DeadChannelMap::IsChannelDead(uint32_t detid, uint32_t channelid, DetectorType dettype, ChannelType chantype) const
    {
        if(!m_isValid)
            return false;

        if (dettype == DetectorType::Barrel1)
        {
            if (chantype == ChannelType::Back)
            {
                channelid += s_sx3BackChannelOffset;
            }
            auto iter = m_barrel1Map.find(GetUUID(detid, channelid));
            if(iter == m_barrel1Map.end())
                return false;
            else
                return true;
        }
        else if (dettype == DetectorType::Barrel2)
        {
            if (chantype == ChannelType::Back)
            {
                channelid += s_sx3BackChannelOffset;
            }
            auto iter = m_barrel2Map.find(GetUUID(detid, channelid));
            if(iter == m_barrel2Map.end())
                return false;
            else
                return true;
        }

        return false;
    }

    bool DeadChannelMap::IsChannelPairDead(uint32_t detid, uint32_t channelfront, uint32_t channelback, DetectorType dettype) const
    {
        if(!m_isValid)
            return false;
            
        if (dettype == DetectorType::Barrel1)
        {
            channelback += s_sx3BackChannelOffset;
            auto iterfront = m_barrel1Map.find(GetUUID(detid, channelfront));
            auto iterback = m_barrel1Map.find(GetUUID(detid, channelback));
            if(iterfront == m_barrel1Map.end() && iterback == m_barrel1Map.end())
                return false;
            else
                return true;
        }
        else if (dettype == DetectorType::Barrel2)
        {
            channelback += s_sx3BackChannelOffset;
            auto iterfront = m_barrel2Map.find(GetUUID(detid, channelfront));
            auto iterback = m_barrel2Map.find(GetUUID(detid, channelback));
            if(iterfront == m_barrel2Map.end() && iterback == m_barrel2Map.end())
                return false;
            else
                return true;
        }

        return false;
    }
}