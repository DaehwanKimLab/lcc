# Copyright 2015 gRPC authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""The Python implementation of the GRPC helloworld.Greeter client."""

from __future__ import print_function

import asyncio
import logging

import sys
sys.path.insert(0, './protos')

import grpc
import lccsimulation_pb2
import lccsimulation_pb2_grpc

import time

async def sim_run() -> None:
    async with grpc.aio.insecure_channel('localhost:50051') as channel:
        stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
        response_iterator = stub.Run(lccsimulation_pb2.Empty())

        prevt = time.time()
        elapsed_time = 0
        count = 0
        async for response in response_iterator:
            curt = time.time()
            deltat = curt - prevt

            prevt = curt
            elapsed_time += deltat

            if elapsed_time >= 1:
                
                print("Received %d sim data in 1 second" % (count))

                elapsed_time = 0
                count = 0

            #for key, val in response.State.Objects.items():
                #print(val.ID, val.Position)
            count += 1

def run():
    while(True):
        cmd = input("Server Function To Run: ")
        if cmd == "Initialize" or cmd == "0":
            with grpc.insecure_channel('localhost:50051') as channel:
                stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
                response = stub.Initialize(lccsimulation_pb2.Empty())
                print(response)

        elif cmd == "Run" or cmd == "1":
            asyncio.run(sim_run())

        elif cmd == "Pause" or cmd == "2":
            with grpc.insecure_channel('localhost:50051') as channel:
                stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
                response = stub.Pause(lccsimulation_pb2.Empty())

        elif cmd == "Stop" or cmd == "3":
            with grpc.insecure_channel('localhost:50051') as channel:
                stub = lccsimulation_pb2_grpc.LCCSimulationStub(channel)
                response = stub.Stop(lccsimulation_pb2.Empty())
        
        elif cmd == "Quit" or cmd == "-1":
            sys.exit(0)

        else:
            print("Function Not Found.")


if __name__ == '__main__':
    logging.basicConfig()
    run()
