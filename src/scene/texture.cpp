
#include "texture.h"

#include <iostream>

namespace Textures {


Spectrum sample_nearest(HDR_Image const &image, Vec2 uv) {
	//clamp texture coordinates, convert to [0,w]x[0,h] pixel space:
	float x = image.w * std::clamp(uv.x, 0.0f, 1.0f);
	float y = image.h * std::clamp(uv.y, 0.0f, 1.0f);

	//the pixel with the nearest center is the pixel that contains (x,y):
	int32_t ix = int32_t(std::floor(x));
	int32_t iy = int32_t(std::floor(y));

	//texture coordinates of (1,1) map to (w,h), and need to be reduced:
	ix = std::min(ix, int32_t(image.w) - 1);
	iy = std::min(iy, int32_t(image.h) - 1);

	return image.at(ix, iy);
}

Spectrum sample_bilinear(HDR_Image const &image, Vec2 uv) {
	// A1T6: sample_bilinear
	//TODO: implement bilinear sampling strategy on texture 'image'
    float x = image.w * std::clamp(uv.x, 0.0f, 1.0f) - 0.5f;
    float y = image.h * std::clamp(uv.y, 0.0f, 1.0f) - 0.5f;
//    printf("\nthe size of image is %d %d", image.w, image.h);
//    printf("\nin sample bilinear uv.x, uv.y are %f %f", uv.x, uv.y);
//    printf("\nin sample bilinear x, y are %f %f", x, y);

    // Find the coordinates of the top-left texel:
    int32_t x0 = int32_t(std::floor(x));
    int32_t y0 = int32_t(std::floor(y));

    // Ensure the coordinates are within the valid range because need x+1, y+1:
    x0 = std::max(0, std::min(x0, int32_t(image.w) - 1));
    y0 = std::max(0, std::min(y0, int32_t(image.h) - 1));

    // Calculate the weights for the interpolation:
    float dx = (float) (x - x0);
    float dy = (float) (y - y0);

    bool isAtVerticalEdge = (uv.x <= 0.5 / image.w || uv.x >= (image.w - 0.5));
    bool isAtHorizontalEdge = (uv.y <= 0.5 / image.h || uv.y >= (image.h - 0.5) / image.h);

    // doesn't work well on right(between R&B)
    if (isAtVerticalEdge && isAtHorizontalEdge) {
        return sample_nearest(image, uv);
    } else if (isAtVerticalEdge || isAtHorizontalEdge) {
        if (isAtVerticalEdge) {
//            printf("\nVertical interpolation");
            Spectrum c0 = image.at(x0, y0);
            Spectrum c1 = y0 < ((int32_t) image.h - 1) ? image.at(x0, y0 + 1) : c0;
            return c0 * (1 - dy) + c1 * dy;
        } else {
//            printf("\nHorizontal interpolation");
            Spectrum c0 = image.at(x0, y0);
            Spectrum c1 = x0 < ((int32_t) image.w - 1) ? image.at(x0 + 1, y0) : c0;
            return c0 * (1 - dx) + c1 * dx;
        }
    }
    else {
        // Get the four texels surrounding the coordinates:
        printf("\nin sample bilinear x0, y0 are %d %d\n", x0, y0);
        Spectrum c00 = image.at(x0, y0);
        Spectrum c10 = x0 < ((int32_t) image.w - 1) ? image.at(x0 + 1, y0) : c00;
        Spectrum c01 = y0 < ((int32_t) image.h - 1) ? image.at(x0, y0 + 1) : c00;
        Spectrum c11 = (x0 < ((int32_t) image.w - 1) && y0 < ((int32_t) image.h - 1)) ? image.at(x0 + 1, y0 + 1) : c00;

        // Perform bilinear interpolation:
        Spectrum cx0 = c00 * (1 - dx) + c10 * dx; // Interpolate between top texels
        Spectrum cx1 = c01 * (1 - dx) + c11 * dx; // Interpolate between bottom texels
        Spectrum cxy = cx0 * (1 - dy) + cx1 * dy; // Interpolate between the results

        return cxy;
    }
}

    Spectrum sample_trilinear(HDR_Image const &base, std::vector< HDR_Image > const &levels, Vec2 uv, float lod) {
	// A1T6: sample_trilinear
	//TODO: implement trilinear sampling strategy on using mip-map 'levels'
//    printf("\n lod is %f", lod);
    int lod_level = std::clamp(int(floor(lod)), 0, int(levels.size()-1));
//    printf("\n lod is %d", lod_level);
    float dlod = lod-float(lod_level);

//    printf("\ncase lower level %d", lod_level-1);
    Spectrum texture_lower = sample_bilinear(lod_level==0 ? base : levels[lod_level-1], uv);
//    printf("\ncase upeer level %d", lod_level);
    Spectrum texture_upper = sample_bilinear(levels[lod_level], uv);
    Spectrum texture = (1-dlod)*texture_lower + dlod*texture_upper;

    return texture;
    // levels[0]会导致Assertion failed: x < w && y < h
//	return sample_bilinear(levels[0], uv); //placeholder so image doesn't look blank
}

/*
 * generate_mipmap- generate mipmap levels from a base image.
 *  base: the base image
 *  levels: pointer to vector of levels to fill (must not be null)
 *
 * generates a stack of levels [1,n] of sizes w_i, h_i, where:
 *   w_i = max(1, floor(w_{i-1})/2)
 *   h_i = max(1, floor(h_{i-1})/2)
 *  with:
 *   w_0 = base.w
 *   h_0 = base.h
 *  and n is the smalles n such that w_n = h_n = 1
 *
 * each level should be calculated by downsampling a blurred version
 * of the previous level to remove high-frequency detail.
 *
 */
void generate_mipmap(HDR_Image const &base, std::vector< HDR_Image > *levels_) {
	assert(levels_);
	auto &levels = *levels_;


	{ // allocate sublevels sufficient to scale base image all the way to 1x1:
		int32_t num_levels = static_cast<int32_t>(std::log2(std::max(base.w, base.h)));
		assert(num_levels >= 0);

		levels.clear();
		levels.reserve(num_levels);

		uint32_t width = base.w;
		uint32_t height = base.h;
		for (int32_t i = 0; i < num_levels; ++i) {
			assert(!(width == 1 && height == 1)); //would have stopped before this if num_levels was computed correctly

			width = std::max(1u, width / 2u);
			height = std::max(1u, height / 2u);

			levels.emplace_back(width, height);
		}
		assert(width == 1 && height == 1);
		assert(levels.size() == uint32_t(num_levels));
	}

	//now fill in the levels using a helper:
	//downsample:
	// fill in dst to represent the low-frequency component of src
        auto downsample = [](HDR_Image const &src, HDR_Image &dst) {
            // dst is half the size of src in each dimension:
            assert(std::max(1u, src.w / 2u) == dst.w);
            assert(std::max(1u, src.h / 2u) == dst.h);

            // A1T6: generate
            // TODO: Write code to fill the levels of the mipmap hierarchy by downsampling

            // Loop over each pixel in dst and compute its value
            for (uint32_t y = 0; y < dst.h; ++y) {
                for (uint32_t x = 0; x < dst.w; ++x) {
                    // Map the pixel in dst back to a 2x2 block in src
                    uint32_t src_x = x * 2;
                    uint32_t src_y = y * 2;

                    Spectrum sum = Spectrum(0.0f);
                    int count = 0;

                    // Sum the color values of the four surrounding pixels in src
                    for (uint32_t dy = 0; dy < 2 && (src_y + dy) < src.h; ++dy) {
                        for (uint32_t dx = 0; dx < 2 && (src_x + dx) < src.w; ++dx) {
                            sum += src.at(src_x + dx, src_y + dy);
                            count++;
                        }
                    }

                    // Calculate the average color value and assign it to the pixel in dst
                    if (count > 0) {
                        dst.at(x, y) = sum / float(count);
                    }
                }
            }
        };

	std::cout << "Regenerating mipmap (" << levels.size() << " levels): [" << base.w << "x" << base.h << "]";
	std::cout.flush();
	for (uint32_t i = 0; i < levels.size(); ++i) {
		HDR_Image const &src = (i == 0 ? base : levels[i-1]);
		HDR_Image &dst = levels[i];
		std::cout << " -> [" << dst.w << "x" << dst.h << "]"; std::cout.flush();

		downsample(src, dst);
	}
	std::cout << std::endl;
	
}

Image::Image(Sampler sampler_, HDR_Image const &image_) {
	sampler = sampler_;
	image = image_.copy();
	update_mipmap();
}

Spectrum Image::evaluate(Vec2 uv, float lod) const {
	if (image.w == 0 && image.h == 0) return Spectrum();
	if (sampler == Sampler::nearest) {
		return sample_nearest(image, uv);
	} else if (sampler == Sampler::bilinear) {
		return sample_bilinear(image, uv);
	} else {
		return sample_trilinear(image, levels, uv, lod);
	}
}

void Image::update_mipmap() {
	if (sampler == Sampler::trilinear) {
		generate_mipmap(image, &levels);
	} else {
		levels.clear();
	}
}

GL::Tex2D Image::to_gl() const {
	return image.to_gl(1.0f);
}

void Image::make_valid() {
	update_mipmap();
}

Spectrum Constant::evaluate(Vec2 uv, float lod) const {
	return color * scale;
}

} // namespace Textures
bool operator!=(const Textures::Constant& a, const Textures::Constant& b) {
	return a.color != b.color || a.scale != b.scale;
}

bool operator!=(const Textures::Image& a, const Textures::Image& b) {
	return a.image != b.image;
}

bool operator!=(const Texture& a, const Texture& b) {
	if (a.texture.index() != b.texture.index()) return false;
	return std::visit(
		[&](const auto& data) { return data != std::get<std::decay_t<decltype(data)>>(b.texture); },
		a.texture);
}
