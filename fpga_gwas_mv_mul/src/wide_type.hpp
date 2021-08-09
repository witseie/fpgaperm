/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef XF_BLAS_TYPES_HPP
#define XF_BLAS_TYPES_HPP

#include <stdint.h>
#include <ostream>
#include <iomanip>
#include <iostream>

#include "hls_math.h"
#include "hls_stream.h"
#include "ap_int.h"
#include "ap_shift_reg.h"

namespace xf
{
    namespace blas
    {
        template <typename T, unsigned int t_Width, unsigned int t_DataWidth = sizeof(T) * 8>
        class WideType
        {
        private:
            T m_Val[t_Width];
            static const unsigned int t_4k = 4096;
            static const unsigned int FLOAT_WIDTH = 7;

        public:
            static const unsigned int t_TypeWidth = t_Width * t_DataWidth;
            typedef ap_int<t_TypeWidth> t_TypeInt; //Fix this for unsigned and signed
            typedef T DataType;
            static const unsigned int t_WidthS = t_Width;
            static const unsigned int t_per4k = t_4k / t_DataWidth / t_Width * 8;

        public:
            T &getVal(unsigned int i)
            {
#ifndef __SYNTHESIS__
                assert(i < t_Width);
#endif
                return (m_Val[i]);
            }

            T &operator[](unsigned int p_Idx)
            {
#ifndef __SYNTHESIS__
                assert(p_Idx < t_Width);
#endif
                return (m_Val[p_Idx]);
            }

            const T &operator[](unsigned int p_Idx) const
            {
#ifndef __SYNTHESIS__
                assert(p_Idx < t_Width);
#endif
                return (m_Val[p_Idx]);
            }

            T *getValAddr()
            {
                return (&m_Val[0]);
            }

            WideType(){
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable = m_Val complete dim = 1
            }

            WideType(const WideType &wt)
            {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable = m_Val complete dim = 1
                for (int i = 0; i < t_Width; i++)
                {
#pragma HLS UNROLL
                    m_Val[i] = wt[i];
                }
            }

            WideType(const t_TypeInt &p_val)
            {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable = m_Val complete dim = 1
                for (int i = 0; i < t_Width; ++i)
                {
#pragma HLS UNROLL
                    ap_uint<t_DataWidth> l_val = p_val.range(t_DataWidth * (1 + i) - 1, t_DataWidth * i);
                    m_Val[i] = *reinterpret_cast<T *>(&l_val);
                }
            }

            WideType(const T p_initScalar)
            {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable = m_Val complete dim = 1
                for (int i = 0; i < t_Width; ++i)
                {
#pragma HLS UNROLL
                    m_Val[i] = p_initScalar;
                }
            }

            operator const t_TypeInt()
            {
#pragma HLS ARRAY_PARTITION variable = m_Val complete dim = 1
                t_TypeInt l_fVal;
                for (int i = 0; i < t_Width; ++i)
                {
#pragma HLS UNROLL
                    T l_v = m_Val[i];
                    ap_uint<t_DataWidth> l_val = *reinterpret_cast<ap_uint<t_DataWidth> *>(&l_v);
                    l_fVal.range(t_DataWidth * (1 + i) - 1, t_DataWidth * i) = l_val;
                }
                return l_fVal;
            }

            T shift(T p_ValIn)
            {
#pragma HLS INLINE
                T l_valOut = m_Val[t_Width - 1];
            WIDE_TYPE_SHIFT:
                for (int i = t_Width - 1; i > 0; --i)
                {
#pragma HLS UNROLL
                    T l_val = m_Val[i - 1];
                    m_Val[i] = l_val;
                }
                m_Val[0] = p_ValIn;
                return (l_valOut);
            }

            T shift()
            {
#pragma HLS INLINE
                T l_valOut = m_Val[t_Width - 1];
                for (int i = t_Width - 1; i > 0; --i)
                {
#pragma HLS UNROLL
                    T l_val = m_Val[i - 1];
                    m_Val[i] = l_val;
                }
                return (l_valOut);
            }

            T unshift()
            {
#pragma HLS INLINE
                T l_valOut = m_Val[0];
                for (int i = 0; i < t_Width - 1; ++i)
                {
#pragma HLS UNROLL
                    T l_val = m_Val[i + 1];
                    m_Val[i] = l_val;
                }
                return (l_valOut);
            }

            T unshift(const T p_val)
            {
#pragma HLS INLINE
                T l_valOut = m_Val[0];
                for (int i = 0; i < t_Width - 1; ++i)
                {
#pragma HLS UNROLL
                    T l_val = m_Val[i + 1];
                    m_Val[i] = l_val;
                }
                m_Val[t_Width - 1] = p_val;
                return (l_valOut);
            }

            static const WideType zero()
            {
                WideType l_zero;
                for (int i = 0; i < t_Width; ++i)
                {
#pragma HLS UNROLL
                    l_zero[i] = 0;
                }
                return (l_zero);
            }

            static unsigned int per4k()
            {
                return (t_per4k);
            }

            void print(std::ostream &os)
            {
                for (int i = 0; i < t_Width; ++i)
                {
                    os << std::setw(FLOAT_WIDTH) << m_Val[i] << " ";
                }
            }

            friend std::ostream &operator<<(std::ostream &os, WideType &p_Val)
            {
                p_Val.print(os);
                return (os);
            }
        };
    }
}

#endif
