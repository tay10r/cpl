#pragma once

#include <memory>

template<typename Texel>
class texture final
{
public:
  texture() = default;

  texture(const int w, const int h, const Texel init)
    : m_data(new Texel[w * h]{ init })
    , m_width(w)
    , m_height(h)
  {
  }

  Texel& operator()(const int x, const int y) { return m_data.get()[(y * m_width) + x]; }

  const Texel& operator()(const int x, const int y) const { return m_data.get()[(y * m_width) + x]; }

  Texel& operator[](const int i) { return m_data.get()[i]; }

  const Texel& operator[](const int i) const { return m_data.get()[i]; }

  int width() const { return m_width; }

  int height() const { return m_height; }

private:
  std::unique_ptr<Texel[]> m_data;

  int m_width{ 0 };

  int m_height{ 0 };
};
