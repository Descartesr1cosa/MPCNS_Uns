// mpcns_pre <setup-file | -> [name=value ...]
#include <cstdio>
#include <exception>
#include <string>

#include "Pipeline.h"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::printf("usage: %s <setup-file | -> [name=value ...]\n"
                    "  task=generate|inspect|export  (see Preprocess/INPUT/setup.txt)\n", argv[0]);
        return 2;
    }
    try
    {
        uns::Params p;
        if (std::string(argv[1]) != "-") p.load_file(argv[1]);
        for (int i = 2; i < argc; ++i) p.set_override(argv[i]);
        const std::string task = p.get_string("task", "generate");
        if (task == "generate") return pre::run_generate(p);
        if (task == "inspect") return pre::run_inspect(p);
        if (task == "export") return pre::run_export(p);
        std::fprintf(stderr, "unknown task '%s'\n", task.c_str());
        return 2;
    }
    catch (const std::exception &e)
    {
        std::fprintf(stderr, "mpcns_pre: error: %s\n", e.what());
        return 1;
    }
}
