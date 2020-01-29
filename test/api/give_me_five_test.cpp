#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "give_me_five.hpp"

using namespace ::testing;

TEST(group1, is_it_five)
{
    int number = my_app::give_me_five();
    EXPECT_EQ(number, 5);
}

TEST(group1, is_it_not_four)
{
    int number = my_app::give_me_five();
    EXPECT_NE(number, 4);
}

